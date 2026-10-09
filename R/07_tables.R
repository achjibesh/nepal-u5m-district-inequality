## ---------------------------------------------------------------------------
## 07_tables.R
##
## Assembles every display table as a tidy, publication-ready data frame and
## writes it to results/tables/ as CSV. manuscript.qmd renders these directly,
## so the numbers in the text and the numbers in the tables cannot drift apart.
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))

DIR_TAB <- file.path(DIR_RESULTS, "tables")
dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)
wr <- function(x, name) {
  write.csv(x, file.path(DIR_TAB, paste0(name, ".csv")), row.names = FALSE)
  invisible(x)
}

births   <- readRDS(file.path(DIR_DERIVED, "births.rds"))
pp       <- readRDS(file.path(DIR_DERIVED, "person_period.rds"))
strat    <- read.csv(file.path(DIR_RESULTS, "rates_by_stratifier.csv"))
nat_per  <- read.csv(file.path(DIR_RESULTS, "rates_national_period.csv"))
valid    <- read.csv(file.path(DIR_RESULTS, "validation_national.csv"))
moran    <- read.csv(file.path(DIR_RESULTS, "moran_global.csv"))
scan_c   <- read.csv(file.path(DIR_RESULTS, "scan_clusters.csv"))
getis    <- read.csv(file.path(DIR_RESULTS, "getis_gi.csv"))
perf     <- read.csv(file.path(DIR_RESULTS, "ml_performance.csv"))
perf_ch  <- read.csv(file.path(DIR_RESULTS, "ml_perf_child.csv"))
spat_cv  <- read.csv(file.path(DIR_RESULTS, "ml_spatial_cv.csv"))
shap_imp <- read.csv(file.path(DIR_RESULTS, "shap_importance.csv"))
inla_cmp <- read.csv(file.path(DIR_RESULTS, "inla_model_comparison.csv"))
inla_fix <- read.csv(file.path(DIR_RESULTS, "inla_fixed_effects.csv"))
inla_var <- read.csv(file.path(DIR_RESULTS, "inla_variance.csv"))
inla_ind <- read.csv(file.path(DIR_RESULTS, "inla_individual_hyper.csv"))
conv     <- read.csv(file.path(DIR_RESULTS, "inla_convergence.csv"))
rr_tab   <- read.csv(file.path(DIR_RESULTS, "inla_district_rr.csv"))
flow     <- read.csv(file.path(DIR_RESULTS, "sample_flow.csv"))

fmt <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")

## Shapefile names carry the post-2015 split as a suffix ("Rukum_e"); spell it
## out for display.
tidy_district <- function(x) {
  x <- gsub("_[eE]\\b", " East", x)
  gsub("_[wW]\\b", " West", x)
}

## ===========================================================================
## Table 1. Sample characteristics and unadjusted under-five mortality
## ===========================================================================

VAR_LABEL <- c(
  sex = "Sex of child", multiple = "Type of birth", birth_order = "Birth order",
  bi_cat = "Preceding birth interval", mage_cat = "Maternal age at birth (years)",
  medu = "Maternal education", ethnicity = "Caste/ethnicity",
  mwork = "Maternal employment",
  media = "Maternal exposure to mass media", wealth = "Household wealth quintile",
  residence = "Place of residence", province = "Province",
  water_imp = "Drinking-water source", sanit_imp = "Sanitation facility",
  electricity = "Household electricity", clean_fuel = "Cooking fuel",
  hh_head_f = "Sex of household head"
)
VAR_ORDER <- names(VAR_LABEL)

dist_tbl <- map_dfr(VAR_ORDER, function(v) {
  births %>%
    filter(!is.na(.data[[v]])) %>%
    group_by(level = as.character(.data[[v]])) %>%
    summarise(n = n(), wn = sum(wt), .groups = "drop") %>%
    mutate(variable = v, pct = 100 * wn / sum(wn)) %>%
    select(variable, level, n, pct)
})

table1 <- dist_tbl %>%
  left_join(strat %>% select(variable, level, u5mr, deaths), by = c("variable", "level")) %>%
  mutate(variable_label = VAR_LABEL[variable],
         variable = factor(variable, levels = VAR_ORDER)) %>%
  arrange(variable) %>%
  ## short column headings: long ones force Word to squeeze the number columns
  transmute(Characteristic = variable_label, Category = level,
            Births = format(n, big.mark = ","),
            `Wt %` = fmt(pct, 1),
            Deaths = deaths,
            U5MR = fmt(u5mr, 1))
wr(table1, "table1_characteristics")

## ===========================================================================
## Table 2. Mortality by period, and validation against the published NDHS
## ===========================================================================

table2a <- nat_per %>%
  transmute(`Birth period` = period,
            `Live births` = format(births, big.mark = ","),
            `Under-five deaths` = deaths,
            `NNMR` = fmt(nnmr, 1), `IMR` = fmt(imr, 1),
            `U5MR (95% CI)` = sprintf("%s (%s-%s)", fmt(u5mr, 1), fmt(lo, 1), fmt(hi, 1)))
wr(table2a, "table2_period_rates")

table2b <- valid %>%
  transmute(`Survey round` = survey_year,
            `Births in the 5 years before the survey` = format(births, big.mark = ","),
            `U5MR, this study` = fmt(u5mr, 1),
            `U5MR, published NDHS report` = published_u5mr,
            `IMR, this study` = fmt(imr, 1),
            `IMR, published` = published_imr,
            `NNMR, this study` = fmt(nnmr, 1),
            `NNMR, published` = published_nnmr)
wr(table2b, "tableS_validation")

## ===========================================================================
## Table 3. Global spatial autocorrelation
## ===========================================================================

table3 <- moran %>%
  filter(rate == "Empirical-Bayes") %>%
  transmute(Stratum = stratum,
            `Moran's I` = sprintf("%.3f", moran_I),
            `p (Moran)` = sprintf("%.4f", moran_p),
            `Geary's C` = sprintf("%.3f", geary_C),
            `p (Geary)` = sprintf("%.4f", geary_p),
            `EB index` = ifelse(is.na(eb_index), "-", sprintf("%.3f", eb_index)),
            `p (EB index)` = ifelse(is.na(eb_p), "-", sprintf("%.4f", eb_p)))
wr(table3, "table3_global_autocorrelation")

## Counts at the nominal threshold and after Benjamini-Hochberg control of
## the 77 simultaneous tests in each stratum.
hot_counts <- bind_rows(
  getis %>% count(stratum, hotspot) %>% mutate(threshold = "p < 0.05"),
  getis %>% count(stratum, hotspot = hotspot_fdr) %>% mutate(threshold = "FDR < 0.05")
) %>%
  filter(hotspot != "Not significant") %>%
  tidyr::pivot_wider(names_from = hotspot, values_from = n, values_fill = 0) %>%
  select(stratum, threshold, everything()) %>%
  arrange(stratum, desc(threshold))
wr(hot_counts, "tableS_getis_counts")

## ===========================================================================
## Table 4. Space-time scan clusters
## ===========================================================================

table4 <- scan_c %>%
  filter(p_value < 0.05 | rank <= 3) %>%
  transmute(Cluster = paste0(rank, ". ", type),
            Period = ifelse(period_from == period_to, period_from,
                            paste0(period_from, " to ", period_to)),
            Districts = n_districts,
            Provinces = provinces,
            Observed = observed, Expected = fmt(expected, 1),
            RR = fmt(rr, 2),
            LLR = fmt(llr, 2),
            ## 999 Monte Carlo replicates: p cannot fall below 1/1000
            `p` = ifelse(p_value <= 0.001, "≤0.001", sprintf("%.3f", p_value)),
            `Member districts` = tidy_district(districts))
wr(table4, "table4_scan_clusters")

## ===========================================================================
## Table 5. Machine-learning performance
## ===========================================================================

table5 <- perf %>%
  transmute(Model = learner,
            `AUROC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", roc_auc, roc_lo, roc_hi),
            `AUPRC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", pr_auc, pr_lo, pr_hi),
            `Brier score` = sprintf("%.5f", brier)) %>%
  left_join(perf_ch %>%
              transmute(Model = learner,
                        `Child-level AUROC` = sprintf("%.3f", roc_auc),
                        `Child-level AUPRC` = sprintf("%.3f", pr_auc)),
            by = "Model") %>%
  mutate(across(everything(), ~ tidyr::replace_na(as.character(.x), "-")))
wr(table5, "table5_ml_performance")

table5b <- spat_cv %>%
  transmute(Province = province,
            `Person-periods` = format(n_person_periods, big.mark = ","),
            `Deaths` = deaths,
            `AUROC, gradient boosting` = sprintf("%.3f", roc_xgb),
            `AUPRC, gradient boosting` = sprintf("%.3f", pr_xgb),
            `AUROC, elastic net` = sprintf("%.3f", roc_enet))
wr(table5b, "tableS_spatial_cv")

## ===========================================================================
## Table 6. Determinants: SHAP importance beside the adjusted odds ratios
## ===========================================================================

## The Bayesian coefficient names are mapped back to a readable label.
tidy_term <- function(x) {
  x %>%
    gsub("^scale\\((.*)\\)$", "\\1", .) %>%
    gsub("seg_f", "Age segment: ", .) %>%
    gsub("^sex", "Sex: ", .) %>%
    gsub("^multiple", "Birth type: ", .) %>%
    gsub("^birth_order", "Birth order: ", .) %>%
    gsub("^bi_cat", "Preceding birth interval: ", .) %>%
    gsub("^mage_cat", "Maternal age at birth: ", .) %>%
    gsub("^medu", "Maternal education: ", .) %>%
    gsub("^wealth", "Wealth quintile: ", .) %>%
    gsub("^residence", "Residence: ", .) %>%
    gsub("^period", "Birth period: ", .) %>%
    gsub("^water_imp", "Drinking water: ", .) %>%
    gsub("^sanit_imp", "Sanitation: ", .) %>%
    gsub("^electricity", "Electricity: ", .) %>%
    gsub("^clean_fuel", "Cooking fuel: ", .) %>%
    gsub("^media", "Media exposure: ", .) %>%
    gsub("^mwork", "Maternal employment: ", .) %>%
    gsub("^hh_head_f", "Household head: ", .) %>%
    gsub("^hh_env", "Household environment: ", .) %>%
    gsub("^ethnicity", "Caste/ethnicity: ", .) %>%
    gsub("alt_km", "Altitude (per SD)", .) %>%
    gsub("^parity$", "Children ever born (per SD)", .) %>%
    gsub("hhsize", "Household size (per SD)", .)
}

table6 <- inla_fix %>%
  filter(term != "(Intercept)") %>%
  mutate(Term = tidy_term(term)) %>%
  transmute(Term,
            `Adjusted OR (95% CrI)` = sprintf("%.2f (%.2f-%.2f)", or, or_lo, or_hi),
            `Posterior P(OR > 1)` = sprintf("%.3f", prob_gt_1))
wr(table6, "table6_adjusted_or")

table6b <- shap_imp %>%
  transmute(Predictor = label, Block = block,
            `Mean |SHAP|` = sprintf("%.4f", mean_abs_shap),
            `Share of total (%)` = sprintf("%.1f", 100 * rel))
wr(table6b, "table6b_shap_importance")

## ===========================================================================
## Table 7. Bayesian model comparison, variance components, convergence
## ===========================================================================

table7a <- inla_cmp %>%
  transmute(Model = model, DIC = fmt(dic, 1), `Effective parameters` = fmt(p_eff, 1),
            WAIC = fmt(waic, 1), `Mean log score` = sprintf("%.4f", log_score))
wr(table7a, "table7_model_comparison")

table7b <- inla_var %>%
  transmute(Component = component,
            `Posterior mean` = sprintf("%.3f", mean),
            `95% CrI` = sprintf("%.3f-%.3f", lower, upper))
wr(table7b, "table7b_variance_components")

table7c <- conv %>%
  transmute(`Birth period` = period,
            `National U5MR` = fmt(u5mr_national, 1),
            `SD of log relative risk` = sprintf("%.3f", sd_log_rr),
            `90th/10th percentile ratio` = sprintf("%.2f", ratio_90_10),
            `Highest/lowest district` = sprintf("%.2f", max_min),
            `Districts with P(RR > 1.2) > 0.8` = n_above_1_2,
            `Absolute 90-10 gap (deaths per 1,000)` = fmt(absolute_gap, 1))
wr(table7c, "table7c_convergence")

conv_post <- read.csv(file.path(DIR_RESULTS, "inla_convergence_posterior.csv"))
## signed to 3 decimals, but a value that rounds to zero is printed unsigned
sg3 <- function(x) { s <- sprintf("%+.3f", x); s[as.numeric(s) == 0] <- "0.000"; s }
table7d <- conv_post %>%
  mutate(first = row_number() == 1) %>%
  transmute(`Birth period` = period,
            `Between-district SD of log RR (95% CrI)` =
              sprintf("%.3f (%.3f-%.3f)", sd_log_rr_median, sd_log_rr_lo, sd_log_rr_hi),
            `Change since 2002-2006 (95% CrI)` =
              ifelse(first, "reference",
                     sprintf("%s (%s to %s)", sg3(change_median), sg3(change_lo), sg3(change_hi))),
            `P(narrower than 2002-2006)` = ifelse(first, "-", sprintf("%.2f", p_narrower)))
wr(table7d, "table7d_convergence_posterior")

sens <- read.csv(file.path(DIR_RESULTS, "inla_sensitivity.csv"))
tableS_sens <- sens %>%
  filter(term != "(Intercept)", !grepl("^seg_f", term)) %>%
  mutate(Term = tidy_term(term),
         est = sprintf("%.2f (%.2f-%.2f)", or, or_lo, or_hi)) %>%
  select(Term, analysis, est) %>%
  tidyr::pivot_wider(names_from = analysis, values_from = est) %>%
  mutate(across(everything(), ~ tidyr::replace_na(.x, "-")))
wr(tableS_sens, "tableS_sensitivity")

mor_tbl <- read.csv(file.path(DIR_RESULTS, "inla_mor.csv"))
wr(mor_tbl %>% transmute(Level = level,
                         `Median odds ratio (95% CrI)` =
                           sprintf("%.2f (%.2f-%.2f)", median_odds_ratio, lo, hi)),
   "tableS_mor")

wr(inla_ind %>%
     transmute(Component = component,
               `Posterior mean` = sprintf("%.3f", mean),
               `95% CrI` = sprintf("%.3f-%.3f", `X0.025quant`, `X0.975quant`)),
   "tableS_individual_hyper")

## ===========================================================================
## Table 8. Districts with the highest posterior burden in the final period
## ===========================================================================

last_p <- tail(PERIOD_LABELS, 1)
table8 <- rr_tab %>%
  filter(period == last_p) %>%
  arrange(desc(rr)) %>%
  slice_head(n = 15) %>%
  transmute(District = tidy_district(district), Province = province,
            `Sampled births` = births, `Observed deaths` = deaths,
            `Raw U5MR` = ifelse(is.na(u5mr_raw), "-", fmt(u5mr_raw, 1)),
            `Posterior U5MR (95% CrI)` = sprintf("%s (%s-%s)",
                                                 fmt(u5mr_smooth, 1),
                                                 fmt(u5mr_smooth_lo, 1),
                                                 fmt(u5mr_smooth_hi, 1)),
            `Relative risk` = fmt(rr, 2),
            `P(RR > 1.2)` = sprintf("%.2f", p_exceed_1_2))
wr(table8, "table8_highest_burden")

## Sample flow --------------------------------------------------------------
wr(flow %>% transmute(Step = step, `Records (n)` = format(n, big.mark = ",")),
   "tableS_sample_flow")

## Key numbers for the narrative -------------------------------------------
## Everything quoted in the Abstract and Results is written here so the text
## can be regenerated mechanically.
key <- list(
  n_births        = nrow(births),
  n_deaths        = sum(births$u5_death, na.rm = TRUE),
  n_person_period = nrow(pp),
  n_clusters      = length(unique(paste(births$survey_year, births$v001))),
  n_districts     = length(unique(births$dcode)),
  u5mr_first      = nat_per$u5mr[1],
  u5mr_last       = nat_per$u5mr[nrow(nat_per)],
  pct_decline     = 100 * (nat_per$u5mr[1] - nat_per$u5mr[nrow(nat_per)]) / nat_per$u5mr[1],
  moran_pooled    = moran$moran_I[moran$stratum == "2002-2022 (pooled)" &
                                    moran$rate == "Empirical-Bayes"],
  moran_pooled_p  = moran$moran_p[moran$stratum == "2002-2022 (pooled)" &
                                    moran$rate == "Empirical-Bayes"],
  scan_rr         = scan_c$rr[1], scan_p = scan_c$p_value[1],
  scan_districts  = scan_c$districts[1],
  best_inla       = inla_cmp$model[1],
  best_auc        = max(perf$roc_auc),
  best_model      = perf$learner[which.max(perf$roc_auc)],
  ## leading predictors excluding the baseline age segment, which necessarily
  ## dominates a hazard model and is not a determinant
  top_shap        = paste(head(shap_imp$label[shap_imp$variable != "seg_f"], 5),
                          collapse = "; ")
)
wr(tibble::tibble(quantity = names(key),
                  value = vapply(key, function(z) as.character(z)[1], character(1))),
   "key_numbers")

message("[07] tables written to ", DIR_TAB)
