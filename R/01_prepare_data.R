## ---------------------------------------------------------------------------
## 01_prepare_data.R
##
## Builds the child-level analytic file and the person-period (discrete-time
## survival) file from the NDHS 2011 / 2016 / 2022 birth recodes, and attaches
## each DHS cluster to one of Nepal's 77 districts by spatial join.
##
## Outputs (data_derived/):
##   births.rds        one row per live birth in the 120 months before interview
##   person_period.rds one row per child x age-segment at risk
##   clusters.rds      cluster -> district / province / altitude lookup
##   district_sf.rds   77-district polygon layer
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))

## ===========================================================================
## 1. Geography: district polygons and DHS cluster locations
## ===========================================================================

district_sf <- st_read(SHP_DISTRICT, quiet = TRUE) |>
  st_transform(4326) |>
  st_make_valid() |>
  transmute(
    dcode    = as.integer(dcode),
    district = str_to_title(str_trim(DISTRICT)),
    prov_shp = as.integer(FIRST_STAT)
  )
stopifnot(nrow(district_sf) == 77)

province_sf <- st_read(SHP_PROVINCE, quiet = TRUE) |>
  st_transform(4326) |> st_make_valid()

## Province names (2015 federal structure, NDHS 2022 nomenclature)
PROV_NAMES <- c("Koshi", "Madhesh", "Bagmati", "Gandaki",
                "Lumbini", "Karnali", "Sudurpashchim")

read_gps <- function(path, round) {
  g <- st_read(path, quiet = TRUE) |> st_transform(4326)
  g |>
    transmute(
      survey_year = as.integer(round),
      v001        = as.integer(DHSCLUST),
      urban       = ifelse(URBAN_RURA == "U", 1L, 0L),
      alt_dem     = suppressWarnings(as.numeric(ALT_DEM)),
      adm1        = as.character(ADM1NAME),
      lat         = as.numeric(LATNUM),
      lon         = as.numeric(LONGNUM)
    ) |>
    filter(!(lat == 0 & lon == 0))          # drop clusters with missing GPS
}

clusters_sf <- map2(GE_FILES, names(GE_FILES), read_gps) |> bind_rows()

## Spatial join to districts. DHS displaces cluster coordinates (up to 2 km
## urban / 5 km rural, 10 km for 1% of rural clusters), so points falling
## marginally outside every polygon are assigned to the nearest district.
joined <- st_join(clusters_sf, district_sf, join = st_within)
miss   <- which(is.na(joined$dcode))
if (length(miss) > 0) {
  nearest <- st_nearest_feature(joined[miss, ], district_sf)
  joined$dcode[miss]    <- district_sf$dcode[nearest]
  joined$district[miss] <- district_sf$district[nearest]
  joined$prov_shp[miss] <- district_sf$prov_shp[nearest]
}
message("[01] clusters joined: ", nrow(joined),
        " (nearest-feature fallback used for ", length(miss), ")")

## Cluster -> district / province lookup (geometry no longer needed). Every
## geographic attribute of a child comes from its district; ecological belt is
## not used in this study.
clusters <- joined |>
  st_drop_geometry() |>
  mutate(province = factor(PROV_NAMES[prov_shp], levels = PROV_NAMES))

district_sf <- district_sf |>
  mutate(province = factor(PROV_NAMES[prov_shp], levels = PROV_NAMES))
stopifnot(!anyNA(clusters$province), !anyNA(district_sf$province))

saveRDS(clusters,    file.path(DIR_DERIVED, "clusters.rds"))
saveRDS(district_sf, file.path(DIR_DERIVED, "district_sf.rds"))
saveRDS(province_sf, file.path(DIR_DERIVED, "province_sf.rds"))

## ===========================================================================
## 2. Birth histories
## ===========================================================================

BR_VARS <- c(
  ## identifiers and sample design
  "caseid", "midx", "v001", "v002", "v005", "v008", "v021", "v022", "v023", "v024",
  ## child
  "b0", "b1", "b2", "b3", "b4", "b5", "b7", "b11", "bord",
  ## mother
  "v011", "v012", "v106", "v130", "v131", "v133", "v201", "v212", "v213", "v501",
  "v502", "v714",
  "v157", "v158", "v159", "v190", "v191",
  ## household and environment
  "v025", "v113", "v116", "v119", "v127", "v136", "v151", "v161", "v745a",
  ## maternal and newborn care (recent births only, used descriptively)
  "m14", "m15", "m17", "m18", "m19", "m3a", "m45", "m70"
)

safe_select <- function(df, vars) {
  present <- intersect(vars, names(df))
  absent  <- setdiff(vars, names(df))
  if (length(absent)) message("   absent in this round: ", paste(absent, collapse = ", "))
  out <- df[, present, drop = FALSE]
  for (v in absent) out[[v]] <- NA
  out[, vars, drop = FALSE]
}

read_br <- function(path, round) {
  message("[01] reading BR ", round)
  d <- read_dta(path)
  d <- safe_select(d, BR_VARS)
  d <- as.data.frame(lapply(d, function(x) {
    if (inherits(x, "haven_labelled")) haven::zap_labels(x) else x
  }), stringsAsFactors = FALSE)
  d$survey_year <- as.integer(round)
  d
}

br_raw <- map2(BR_FILES, names(BR_FILES), read_br) |> bind_rows()
message("[01] pooled birth records: ", nrow(br_raw))

## ---------------------------------------------------------------------------
## 2.1 Child-level analytic file
## ---------------------------------------------------------------------------

births <- br_raw |>
  mutate(
    age_at_int = v008 - b3,                 # completed months birth -> interview
    birth_ad   = bs_cmc_to_ad(b3),          # Gregorian decimal year
    birth_year = floor(birth_ad),
    int_ad     = bs_cmc_to_ad(v008),
    wt         = v005 / 1e6,
    psu        = v021,
    strata     = if_else(is.na(v023), v022, v023)
  ) |>
  filter(age_at_int >= 0, age_at_int < WINDOW_MONTHS) |>
  mutate(period = cut(birth_year, breaks = PERIOD_BREAKS, right = FALSE,
                      labels = PERIOD_LABELS)) |>
  filter(!is.na(period))

## Outcome and censoring.
## b5 = 1 alive / 0 dead; b7 = age at death in completed months.
births <- births |>
  mutate(
    dead      = if_else(b5 == 0, 1L, 0L),
    age_death = if_else(dead == 1L, as.numeric(b7), NA_real_),
    
    ## true exit time from the 0-59 month window (months since birth)
    exit_month = case_when(
      dead == 1L & !is.na(age_death) ~ pmin(age_death + 0.0, 60),
      dead == 1L &  is.na(age_death) ~ NA_real_,
      TRUE                           ~ pmin(age_at_int, 60)
    ),
    u5_death = case_when(
      dead == 1L & !is.na(age_death) & age_death <  60 ~ 1L,
      dead == 1L & !is.na(age_death) & age_death >= 60 ~ 0L,
      dead == 0L                                       ~ 0L,
      TRUE                                             ~ NA_integer_
    ),
    complete_exposure = age_at_int >= 60 | u5_death %in% 1L,
    neonatal_death = as.integer(dead == 1L & !is.na(age_death) & age_death <  1),
    infant_death   = as.integer(dead == 1L & !is.na(age_death) & age_death < 12)
  ) |>
  filter(!is.na(exit_month), !is.na(u5_death))

## Upper limit of the at-risk set. A child that dies is at risk in the segment
## in which the death occurs, so half a month is added to its exit time; a
## surviving child is at risk up to its attained age.
births <- births |>
  mutate(risk_until = if_else(u5_death == 1L, exit_month + 0.5, exit_month))

## ---------------------------------------------------------------------------
## 2.2 Covariates
## ---------------------------------------------------------------------------

births <- births |>
  left_join(clusters |> select(survey_year, v001, dcode, district, province,
                               alt_dem, urban_gps = urban, lat, lon),
            by = c("survey_year", "v001")) |>
  mutate(
    ## --- child ---------------------------------------------------------
    sex         = factor(if_else(b4 == 1, "Male", "Female"), levels = c("Female", "Male")),
    multiple    = factor(if_else(b0 == 0, "Singleton", "Multiple"),
                         levels = c("Singleton", "Multiple")),
    birth_order = factor(case_when(bord == 1 ~ "1",
                                   bord %in% 2:3 ~ "2-3",
                                   bord %in% 4:5 ~ "4-5",
                                   bord >= 6 ~ "6+"),
                         levels = c("1", "2-3", "4-5", "6+")),
    bord_num    = as.numeric(bord),
    bi_cat      = factor(case_when(is.na(b11) ~ "First birth",
                                   b11 <  24  ~ "<24 months",
                                   b11 <  36  ~ "24-35 months",
                                   TRUE       ~ ">=36 months"),
                         levels = c(">=36 months", "24-35 months",
                                    "<24 months", "First birth")),
    bi_num      = as.numeric(b11),

    ## --- mother --------------------------------------------------------
    mage_birth  = (b3 - v011) / 12,
    mage_cat    = factor(case_when(mage_birth < 20 ~ "<20",
                                   mage_birth < 25 ~ "20-24",
                                   mage_birth < 30 ~ "25-29",
                                   mage_birth < 35 ~ "30-34",
                                   TRUE            ~ "35+"),
                         levels = c("20-24", "<20", "25-29", "30-34", "35+")),
    medu        = factor(case_when(v106 == 0 ~ "None",
                                   v106 == 1 ~ "Primary",
                                   v106 %in% 2:3 ~ "Secondary or higher"),
                         levels = c("Secondary or higher", "Primary", "None")),
    medu_years  = as.numeric(v133),
    parity      = as.numeric(v201),
    mwork       = factor(if_else(v714 == 1, "Employed", "Not employed"),
                         levels = c("Not employed", "Employed")),
    media       = factor(if_else((v157 %in% 1:3) | (v158 %in% 1:3) | (v159 %in% 1:3),
                                 "Yes", "No"), levels = c("No", "Yes")),

    ## --- household -----------------------------------------------------
    wealth      = factor(case_when(v190 == 1 ~ "Poorest", v190 == 2 ~ "Poorer",
                                   v190 == 3 ~ "Middle",  v190 == 4 ~ "Richer",
                                   v190 == 5 ~ "Richest"),
                         levels = c("Richest", "Richer", "Middle", "Poorer", "Poorest")),
    wealth_score = as.numeric(v191) / 1e5,
    residence   = factor(if_else(v025 == 1, "Urban", "Rural"), levels = c("Urban", "Rural")),
    hhsize      = as.numeric(v136),
    hh_head_f   = factor(if_else(v151 == 2, "Female", "Male"), levels = c("Male", "Female")),

    ## Caste/ethnicity (v131) is coded identically in all three rounds. It is
    ## grouped into the seven categories conventionally used in NDHS further
    ## analyses, with Brahmin/Chhetri as the reference.
    ethnicity   = factor(case_when(v131 %in% 1:3 ~ "Brahmin/Chhetri",
                                   v131 == 4     ~ "Terai/Madhesi other caste",
                                   v131 %in% 5:6 ~ "Dalit",
                                   v131 == 7     ~ "Newar",
                                   v131 %in% 8:9 ~ "Janajati",
                                   v131 == 10    ~ "Muslim",
                                   !is.na(v131)  ~ "Other"),
                         levels = c("Brahmin/Chhetri", "Terai/Madhesi other caste",
                                    "Dalit", "Newar", "Janajati", "Muslim", "Other")),

    ## DHS code 97 ("not a de jure resident") marks a mother who was a visitor
    ## in the sampled household, so that household's water, sanitation, fuel
    ## and electricity were never collected. This affects the same 1,627
    ## children in all four variables, spread across every wealth quintile
    ## (211 of them in the richest). Their weighted under-five mortality lies
    ## between that of improved- and unimproved-water households, so there is
    ## no basis for assigning them to either: folding them into "Unimproved" /
    ## "Polluting" would give 5% of children the worst category on no evidence,
    ## and dropping them would discard a group that is not missing at random.
    ## They are carried as an explicit level instead.
    water_imp   = factor(case_when(
                           v113 == 97 ~ "Not recorded",
                           v113 %in% c(10, 11, 12, 13, 14, 21, 31, 41, 51,
                                       61, 62, 71, 72) ~ "Improved",
                           !is.na(v113) ~ "Unimproved"),
                         levels = c("Improved", "Unimproved", "Not recorded")),
    sanit_imp   = factor(case_when(
                           v116 == 97 ~ "Not recorded",
                           v116 %in% c(10, 11, 12, 13, 14, 15, 21, 22, 41) ~ "Improved",
                           !is.na(v116) ~ "Unimproved"),
                         levels = c("Improved", "Unimproved", "Not recorded")),
    electricity = factor(case_when(v119 == 97 | v119 == 7 ~ "Not recorded",
                                   v119 == 1 ~ "Yes",
                                   !is.na(v119) ~ "No"),
                         levels = c("Yes", "No", "Not recorded")),
    clean_fuel  = factor(case_when(v161 == 97 ~ "Not recorded",
                                   v161 %in% 1:5 ~ "Clean",
                                   !is.na(v161) ~ "Polluting"),
                         levels = c("Clean", "Polluting", "Not recorded")),

    ## --- geography and time --------------------------------------------
    altitude   = as.numeric(alt_dem),
    alt_km     = altitude / 1000,
    round      = factor(survey_year, levels = c(2011, 2016, 2022)),
    period     = factor(period, levels = PERIOD_LABELS),
    period_num = as.integer(period)
  )

## STROBE sample-size accounting ---------------------------------------------
n_window <- br_raw |>
  mutate(a = v008 - b3) |>
  filter(a >= 0, a < WINDOW_MONTHS) |>
  nrow()

flow <- tibble::tibble(
  step = c("Live births reported in birth histories, NDHS 2011 + 2016 + 2022",
           "Born within 120 months preceding the interview",
           "Birth year within 2002-2022 and age at death non-missing",
           "Cluster GPS available and linked to one of 77 districts",
           "Child-level analytic sample"),
  n = c(nrow(br_raw), n_window, nrow(births),
        sum(!is.na(births$dcode)), sum(!is.na(births$dcode)))
)

births <- births |> filter(!is.na(dcode))
births$child_id <- seq_len(nrow(births))

message("[01] child-level analytic sample: ", nrow(births),
        " | under-five deaths: ", sum(births$u5_death, na.rm = TRUE))

write.csv(flow, file.path(DIR_RESULTS, "sample_flow.csv"), row.names = FALSE)

## ===========================================================================
## 3. Person-period (discrete-time survival) expansion
## ===========================================================================
## Each child contributes one row per DHS age segment in which it was at risk.
## `died` = 1 in the segment containing the age at death. Exposure is truncated
## at the attained age at interview, which removes the right-censoring bias
## that arises when a binary under-five death indicator is fitted directly to
## a retrospective birth window.

expand_person_period <- function(dat, segments = AGE_SEGMENTS) {
  dt  <- as.data.table(dat)[, .(child_id, exit_month, risk_until, age_death, u5_death)]
  seg <- as.data.table(segments)
  pp  <- CJ(child_id = dt$child_id, seg = seg$seg)
  pp  <- merge(pp, seg, by = "seg")
  pp  <- merge(pp, dt, by = "child_id")
  pp  <- pp[start < risk_until]                       # entered the segment
  pp[, died := as.integer(u5_death == 1L & !is.na(age_death) &
                            age_death >= start & age_death < end)]
  ## months at risk within the segment (>= 0.5 so that a death on the first
  ## day of a segment still contributes exposure)
  pp[, expo  := pmax(pmin(exit_month, end) - start, 0.5)]
  pp[, seg_f := factor(label, levels = segments$label)]
  setorder(pp, child_id, seg)
  pp[, .(child_id, seg, seg_f, start, end, width, expo, died)]
}

pp <- expand_person_period(births)

person_period <- merge(as.data.table(pp),
                       as.data.table(births)[, !c("exit_month"), with = FALSE],
                       by = "child_id", all.x = TRUE) |>
  as_tibble()

message("[01] person-period rows: ", nrow(person_period),
        " | segment deaths: ", sum(person_period$died))

saveRDS(births,        file.path(DIR_DERIVED, "births.rds"))
saveRDS(person_period, file.path(DIR_DERIVED, "person_period.rds"))

message("[01] done.")

