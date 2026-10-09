## ---------------------------------------------------------------------------
## 02_mortality_rates.R
##
## Life-table estimation (by birth cohort) of neonatal, infant and under-five
## mortality from the person-period file:
##   * national estimates for the five years preceding each survey
##     (validation against the published NDHS key indicators)
##   * national and provincial estimates by calendar birth period
##   * district x period counts and rates that feed the spatial and
##     Bayesian spatio-temporal models
##
## Outputs (results/, data_derived/):
##   validation_national.csv, rates_national_period.csv,
##   rates_province_period.csv, district_period.rds / .csv
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))

births <- readRDS(file.path(DIR_DERIVED, "births.rds"))
pp     <- readRDS(file.path(DIR_DERIVED, "person_period.rds"))

## ===========================================================================
## 1. Life-table machinery
## ===========================================================================
## Segment-specific probability of dying, q_i, uses the standard DHS
## denominator: children entering the segment, with those censored part-way
## through it counted as one half. The cumulative probability of dying before
## exact age a is 1 - prod(1 - q_i) over the segments below a.

seg_probs <- function(d, wcol = "wt") {
  w <- d[[wcol]]
  d %>%
    mutate(.w = w,
           ## a child is censored inside the segment if it neither died there
           ## nor survived to the segment's upper boundary
           censored_in = as.integer(died == 0L & expo < width)) %>%
    group_by(seg, seg_f, start, end) %>%
    summarise(
      deaths   = sum(.w * died),
      entrants = sum(.w),
      cens     = sum(.w * censored_in),
      n_raw    = n(),
      d_raw    = sum(died),
      .groups  = "drop"
    ) %>%
    mutate(denom = entrants - 0.5 * cens,
           q     = if_else(denom > 0, deaths / denom, 0))
}

cum_q <- function(sp, upper) {
  s <- sp %>% filter(end <= upper)
  1 - prod(1 - s$q)
}

## Neonatal (< 1 month), infant (< 12 months) and under-five (< 60 months)
## mortality per 1,000 live births.
mort_rates <- function(d, wcol = "wt") {
  sp <- seg_probs(d, wcol)
  tibble::tibble(
    nnmr   = 1000 * cum_q(sp, 1),
    imr    = 1000 * cum_q(sp, 12),
    u5mr   = 1000 * cum_q(sp, 60),
    births = length(unique(d$child_id)),
    deaths = sum(d$died)
  )
}

## Cluster (PSU) bootstrap confidence interval, resampling primary sampling
## units within survey round to respect the complex design.
boot_ci <- function(d, R = 400, wcol = "wt", stat = "u5mr", seed = SEED) {
  set.seed(seed)
  psus <- unique(d[, c("survey_year", "psu")])
  dt   <- as.data.table(d)
  setkey(dt, survey_year, psu)
  out <- numeric(R)
  for (r in seq_len(R)) {
    idx <- psus[sample.int(nrow(psus), nrow(psus), replace = TRUE), ]
    bs  <- dt[as.data.table(idx), on = .(survey_year, psu), allow.cartesian = TRUE]
    out[r] <- mort_rates(bs, wcol)[[stat]]
  }
  stats::quantile(out, c(0.025, 0.975), na.rm = TRUE)
}

## ===========================================================================
## 2. Validation: five years preceding each survey
## ===========================================================================
## The published NDHS key-indicator estimates for the 0-4 years before the
## survey are 54 (2011), 39 (2016) and 33 (2022) deaths per 1,000 live births.

pp5 <- pp %>% filter(age_at_int < 60)

validation <- pp5 %>%
  group_by(survey_year) %>%
  group_modify(~ mort_rates(.x)) %>%
  ungroup() %>%
  mutate(published_u5mr = c(54, 39, 33)[match(survey_year, c(2011, 2016, 2022))],
         published_imr  = c(46, 32, 28)[match(survey_year, c(2011, 2016, 2022))],
         published_nnmr = c(33, 21, 21)[match(survey_year, c(2011, 2016, 2022))])

print(validation)
write.csv(validation, file.path(DIR_RESULTS, "validation_national.csv"), row.names = FALSE)

## ===========================================================================
## 3. National and provincial estimates by calendar birth period
## ===========================================================================

nat_period <- pp %>%
  group_by(period) %>%
  group_modify(~ mort_rates(.x)) %>%
  ungroup()

nat_period_ci <- pp %>%
  group_split(period) %>%
  map_dfr(function(x) {
    ci <- boot_ci(x, R = 300)
    tibble::tibble(period = unique(x$period), lo = ci[1], hi = ci[2])
  })

nat_period <- nat_period %>% left_join(nat_period_ci, by = "period")
print(nat_period)
write.csv(nat_period, file.path(DIR_RESULTS, "rates_national_period.csv"), row.names = FALSE)

prov_period <- pp %>%
  group_by(province, period) %>%
  group_modify(~ mort_rates(.x)) %>%
  ungroup()
write.csv(prov_period, file.path(DIR_RESULTS, "rates_province_period.csv"), row.names = FALSE)

## Rates by every stratifier, for the descriptive table --------------------
strat_vars <- c("sex", "residence", "wealth", "medu", "ethnicity", "province", "birth_order",
                "bi_cat", "mage_cat", "multiple", "water_imp", "sanit_imp",
                "electricity", "clean_fuel", "media", "mwork", "hh_head_f")

rates_by_strat <- map_dfr(strat_vars, function(v) {
  pp %>%
    filter(!is.na(.data[[v]])) %>%
    group_by(level = as.character(.data[[v]])) %>%
    group_modify(~ mort_rates(.x)) %>%
    ungroup() %>%
    mutate(variable = v, .before = 1)
})
write.csv(rates_by_strat, file.path(DIR_RESULTS, "rates_by_stratifier.csv"), row.names = FALSE)

## ===========================================================================
## 4. District x period file for the spatial models
## ===========================================================================
## For every district-period cell we retain
##   * the weighted life-table U5MR,
##   * unweighted counts of births and under-five deaths, and
##   * the number of deaths expected under indirect standardisation by the
##     national age-segment schedule -- the offset used by the Poisson scan
##     statistic and the Bayesian spatio-temporal model.

## Two sets of expected counts are produced.
##   expected_crude   - indirect standardisation by age segment only. A cluster
##                      detected against this baseline reflects the combined
##                      effect of place and calendar time.
##   expected         - indirect standardisation by age segment *within period*.
##                      This removes the national secular decline, so that a
##                      detected cluster is one whose mortality is high
##                      relative to the rest of Nepal at the same date. It is
##                      the offset used for the scan statistic and the
##                      Bayesian spatio-temporal model.

nat_sp <- seg_probs(pp) %>% select(seg, q_nat = q)          # pooled schedule

nat_sp_period <- pp %>%
  group_by(period) %>%
  group_modify(~ seg_probs(.x)) %>%
  ungroup() %>%
  select(period, seg, q_nat_p = q)

expected_deaths <- pp %>%
  left_join(nat_sp, by = "seg") %>%
  left_join(nat_sp_period, by = c("period", "seg")) %>%
  group_by(dcode, district, province, period) %>%
  summarise(
    births         = n_distinct(child_id),
    deaths         = sum(died),
    expected       = sum(q_nat_p * (expo / width)),
    expected_crude = sum(q_nat   * (expo / width)),
    pmonths        = sum(expo),
    .groups        = "drop"
  )

district_rates <- pp %>%
  group_by(dcode, period) %>%
  group_modify(~ mort_rates(.x)) %>%
  ungroup() %>%
  select(dcode, period, u5mr, imr, nnmr)

district_period <- expected_deaths %>%
  left_join(district_rates, by = c("dcode", "period")) %>%
  mutate(
    smr        = if_else(expected > 0, deaths / expected, NA_real_),
    crude_rate = 1000 * deaths / births,
    period_num = as.integer(factor(period, levels = PERIOD_LABELS))
  )

## Complete the 77 x 4 lattice so that the Bayesian model estimates every cell,
## including districts with no sampled births in a given period.
district_sf <- readRDS(file.path(DIR_DERIVED, "district_sf.rds"))
lattice <- tidyr::expand_grid(
  district_sf %>% st_drop_geometry() %>% select(dcode, district, province),
  period = factor(PERIOD_LABELS, levels = PERIOD_LABELS)
)

district_period <- lattice %>%
  left_join(district_period %>% select(-district, -province),
            by = c("dcode", "period")) %>%
  mutate(
    across(c(births, deaths), ~ tidyr::replace_na(.x, 0)),
    period_num = as.integer(period),
    dist_id    = as.integer(factor(dcode, levels = sort(unique(district_sf$dcode))))
  )

message("[02] district-period cells: ", nrow(district_period),
        " | cells with no sampled births: ", sum(district_period$births == 0))

saveRDS(district_period, file.path(DIR_DERIVED, "district_period.rds"))
write.csv(district_period, file.path(DIR_RESULTS, "district_period.csv"), row.names = FALSE)

## District totals pooled over periods (for the cluster maps) ---------------
district_all <- pp %>%
  group_by(dcode) %>%
  group_modify(~ mort_rates(.x)) %>%
  ungroup()
saveRDS(district_all, file.path(DIR_DERIVED, "district_all.rds"))

message("[02] done.")
