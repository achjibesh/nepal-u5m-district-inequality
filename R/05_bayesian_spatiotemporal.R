## ---------------------------------------------------------------------------
## 05_bayesian_spatiotemporal.R
##
## Two Bayesian hierarchical models fitted with integrated nested Laplace
## approximation (INLA).
##
## Model A - small-area spatio-temporal model on district x period death
##   counts. Observed deaths are Poisson with an offset equal to the number of
##   deaths expected under indirect standardisation by age segment within
##   period, so that a relative risk of 1 means "the national level at that
##   date". The linear predictor carries
##     * a BYM2 district effect (a mixture of an intrinsic conditional
##       autoregressive and an independent component, with the mixing
##       parameter phi giving the share of area variance that is spatially
##       structured),
##     * a first-order random-walk period effect, and
##     * a space-time interaction. All four Knorr-Held interaction types are
##       fitted and compared by DIC, WAIC and the log score.
##   Penalised-complexity priors are used throughout.
##
## Model B - individual-level discrete-time survival model on the person-period
##   file, giving adjusted odds ratios for each determinant together with the
##   district-level spatial variance that survives adjustment. This is the
##   classical epidemiological counterpart to the machine-learning analysis
##   and quantifies how much of the spatial inequity in child survival is left
##   unexplained by measured household and maternal characteristics.
##
## Outputs (results/): inla_model_comparison.csv, inla_district_rr.csv,
##   inla_variance.csv, inla_fixed_effects.csv, inla_convergence.csv
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))
suppressPackageStartupMessages(library(INLA))
INLA::inla.setOption(num.threads = "4:1")

district_sf <- readRDS(file.path(DIR_DERIVED, "district_sf.rds")) %>% arrange(dcode)
dp          <- readRDS(file.path(DIR_DERIVED, "district_period.rds")) %>%
  arrange(period, dcode)
wts         <- readRDS(file.path(DIR_DERIVED, "spatial_weights.rds"))
nat_period  <- read.csv(file.path(DIR_RESULTS, "rates_national_period.csv"))

## ===========================================================================
## 1. Graph
## ===========================================================================
adj_file <- file.path(DIR_DERIVED, "nepal_districts.adj")
nb2INLA(adj_file, wts$nb)
g <- inla.read.graph(adj_file)

dp <- dp %>%
  mutate(
    id_space  = as.integer(factor(dcode, levels = district_sf$dcode)),
    id_time   = as.integer(period),
    id_space2 = id_space,
    id_time2  = id_time,
    id_st     = seq_len(n()),
    ## District-period cells with no sampled births have no expected count.
    ## They are kept in the lattice with a token exposure so that the model
    ## still produces a posterior relative risk for them, interpolated from
    ## neighbouring districts and adjacent periods rather than from data.
    E         = pmax(tidyr::replace_na(expected, 0), 1e-4)
  )
stopifnot(!any(is.na(dp$E)), !any(is.na(dp$deaths)))

## PC priors: P(sd > 1) = 0.01 for the random effects, and a uniform prior on
## the BYM2 mixing parameter shifted towards unstructured heterogeneity.
pc_prec <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))
pc_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                phi  = list(prior = "pc", param = c(0.5, 0.5)))

base_f <- deaths ~ 1 +
  f(id_space, model = "bym2", graph = g, scale.model = TRUE,
    constr = TRUE, hyper = pc_bym2) +
  f(id_time, model = "rw1", scale.model = TRUE, constr = TRUE, hyper = pc_prec)

## Knorr-Held interaction types -------------------------------------------
## I   unstructured in space and time
## II  district-specific random walks (structured time, unstructured space)
## III period-specific spatial fields (structured space, unstructured time)
## IV  fully structured Kronecker interaction
forms <- list(
  `No interaction` = base_f,
  `Type I`   = update(base_f, . ~ . + f(id_st, model = "iid", hyper = pc_prec)),
  `Type II`  = update(base_f, . ~ . +
                        f(id_time2, model = "rw1", group = id_space2,
                          control.group = list(model = "iid"),
                          scale.model = TRUE, hyper = pc_prec)),
  `Type III` = update(base_f, . ~ . +
                        f(id_space2, model = "besag", graph = g, group = id_time2,
                          control.group = list(model = "iid"),
                          scale.model = TRUE, hyper = pc_prec)),
  `Type IV`  = update(base_f, . ~ . +
                        f(id_space2, model = "besag", graph = g, group = id_time2,
                          control.group = list(model = "rw1"),
                          scale.model = TRUE, hyper = pc_prec))
)

fit_one <- function(f) {
  inla(f, family = "poisson", data = dp, E = dp$E,
       control.predictor = list(compute = TRUE, link = 1),
       control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE,
                              config = TRUE, return.marginals.predictor = TRUE),
       control.inla = list(strategy = "adaptive", int.strategy = "ccd"),
       verbose = FALSE)
}

message("[05] fitting area-level models ...")
fits <- lapply(names(forms), function(nm) {
  message("   ", nm)
  tryCatch(fit_one(forms[[nm]]), error = function(e) {
    message("     failed: ", conditionMessage(e)); NULL })
})
names(fits) <- names(forms)
fits <- fits[!vapply(fits, is.null, logical(1))]

comparison <- map_dfr(names(fits), function(nm) {
  m <- fits[[nm]]
  tibble::tibble(model = nm,
                 dic = m$dic$dic, p_eff = m$dic$p.eff, waic = m$waic$waic,
                 log_score = -mean(log(pmax(m$cpo$cpo, 1e-12)), na.rm = TRUE),
                 n_failure = sum(m$cpo$failure > 0, na.rm = TRUE))
}) %>% arrange(dic)
print(as.data.frame(comparison), digits = 5)
write.csv(comparison, file.path(DIR_RESULTS, "inla_model_comparison.csv"), row.names = FALSE)

best_name <- comparison$model[1]
best <- fits[[best_name]]
message("[05] best area model: ", best_name)
saveRDS(list(name = best_name, summary_fixed = best$summary.fixed,
             summary_hyper = best$summary.hyperpar),
        file.path(DIR_DERIVED, "inla_best_area.rds"))

## ===========================================================================
## 2. Posterior district-period relative risks and smoothed rates
## ===========================================================================

fitted_rr <- best$summary.fitted.values[seq_len(nrow(dp)), ]
exceed <- vapply(seq_len(nrow(dp)), function(i)
  1 - INLA::inla.pmarginal(1.2, best$marginals.fitted.values[[i]]), numeric(1))
exceed1 <- vapply(seq_len(nrow(dp)), function(i)
  1 - INLA::inla.pmarginal(1.0, best$marginals.fitted.values[[i]]), numeric(1))

nat_lookup <- setNames(nat_period$u5mr, nat_period$period)

district_rr <- dp %>%
  select(dcode, district, province, period, period_num,
         births, deaths, expected, u5mr_raw = u5mr, u5mr_eb) %>%
  mutate(
    rr       = fitted_rr$mean,
    rr_lo    = fitted_rr$`0.025quant`,
    rr_hi    = fitted_rr$`0.975quant`,
    rr_sd    = fitted_rr$sd,
    p_exceed_1_2 = exceed,
    p_exceed_1   = exceed1,
    u5mr_smooth  = rr * nat_lookup[as.character(period)],
    u5mr_smooth_lo = rr_lo * nat_lookup[as.character(period)],
    u5mr_smooth_hi = rr_hi * nat_lookup[as.character(period)]
  )
write.csv(district_rr, file.path(DIR_RESULTS, "inla_district_rr.csv"), row.names = FALSE)
saveRDS(district_rr, file.path(DIR_DERIVED, "inla_district_rr.rds"))

message("[05] districts with P(RR > 1.2) > 0.8 in the final period: ",
        sum(district_rr$p_exceed_1_2[district_rr$period == tail(PERIOD_LABELS, 1)] > 0.8))

## ===========================================================================
## 3. Variance decomposition and convergence
## ===========================================================================

hyper <- best$summary.hyperpar
sd_from_prec <- function(row_name) {
  if (!row_name %in% rownames(hyper)) return(NA_real_)
  1 / sqrt(hyper[row_name, "mean"])
}
variance_tbl <- tibble::tibble(
  component = rownames(hyper),
  mean = hyper[, "mean"], sd = hyper[, "sd"],
  lower = hyper[, "0.025quant"], upper = hyper[, "0.975quant"]
)
write.csv(variance_tbl, file.path(DIR_RESULTS, "inla_variance.csv"), row.names = FALSE)
print(as.data.frame(variance_tbl), digits = 4)

## Spatial inequality over time: the dispersion of posterior district relative
## risks within each period. A narrowing distribution means districts are
## converging on the national trajectory; a stable or widening one means the
## geography of disadvantage is entrenched.
convergence <- district_rr %>%
  group_by(period) %>%
  summarise(
    n_districts = n(),
    sd_log_rr   = sd(log(rr)),
    iqr_rr      = IQR(rr),
    ratio_90_10 = quantile(rr, 0.9) / quantile(rr, 0.1),
    max_min     = max(rr) / min(rr),
    n_above_1_2 = sum(p_exceed_1_2 > 0.8),
    .groups = "drop"
  ) %>%
  left_join(nat_period %>% select(period, u5mr_national = u5mr), by = "period") %>%
  mutate(absolute_gap = (quantile_hi <- NA_real_))
convergence$absolute_gap <- district_rr %>%
  group_by(period) %>%
  summarise(g = (quantile(rr, 0.9) - quantile(rr, 0.1)) *
              nat_lookup[as.character(first(period))], .groups = "drop") %>%
  pull(g)
print(as.data.frame(convergence), digits = 4)
write.csv(convergence, file.path(DIR_RESULTS, "inla_convergence.csv"), row.names = FALSE)

## The statistics above use posterior means, which are shrunk towards the
## national level by an amount that depends on how much information a period
## carries. The latest period has the fewest births and hence the strongest
## shrinkage, which on its own would make districts look more alike and mimic
## convergence. Dispersion is therefore recomputed inside each of 1,000 joint
## posterior draws, giving a posterior distribution -- with a credible
## interval -- for the between-district SD of the log relative risk in each
## period and for its change since the first period.
set.seed(SEED)
NDRAW <- 1000
draws <- INLA::inla.posterior.sample(NDRAW, best)
pred_idx <- grep("^Predictor:", rownames(draws[[1]]$latent))
stopifnot(length(pred_idx) == nrow(dp))
eta <- vapply(draws, function(s) s$latent[pred_idx, 1], numeric(nrow(dp)))
disp_draws <- purrr::map_dfr(PERIOD_LABELS, function(p) {
  rows <- which(dp$period == p)
  tibble::tibble(period = p, draw = seq_len(NDRAW),
                 sd_log_rr = apply(eta[rows, , drop = FALSE], 2, sd))
})
disp_change <- disp_draws %>%
  left_join(disp_draws %>% filter(period == PERIOD_LABELS[1]) %>%
              select(draw, sd_first = sd_log_rr), by = "draw") %>%
  mutate(change = sd_log_rr - sd_first)
convergence_post <- disp_change %>%
  group_by(period) %>%
  summarise(sd_log_rr_median = median(sd_log_rr),
            sd_log_rr_lo = quantile(sd_log_rr, 0.025),
            sd_log_rr_hi = quantile(sd_log_rr, 0.975),
            change_median = median(change),
            change_lo = quantile(change, 0.025),
            change_hi = quantile(change, 0.975),
            p_narrower = mean(change < 0),
            .groups = "drop") %>%
  mutate(period = factor(period, levels = PERIOD_LABELS)) %>%
  arrange(period)
print(as.data.frame(convergence_post), digits = 3)
write.csv(convergence_post, file.path(DIR_RESULTS, "inla_convergence_posterior.csv"),
          row.names = FALSE)

## ===========================================================================
## 3b. Cluster detection on the posterior relative-risk surface
## ===========================================================================
## Exploratory only: these outputs are NOT used as evidence in the manuscript.
## The posterior surface is smoothed by the model's own spatial (BYM2) prior,
## so spatial autocorrelation measured on it is partly created by the model;
## testing for clustering on it would be circular. The evidence for clustering
## comes from script 03 (empirical-Bayes rates, scan statistic) and from the
## BYM2 mixing parameter. The maps are kept as a descriptive check.

lw <- wts$listw
nb_star <- spdep::include.self(wts$nb)
lw_star <- spdep::nb2listw(nb_star, style = "B", zero.policy = TRUE)

lisa_post <- map_dfr(PERIOD_LABELS, function(p) {
  d  <- district_rr %>% filter(period == p) %>% arrange(dcode)
  x  <- log(d$rr)
  lm_res <- spdep::localmoran_perm(x, lw, nsim = 9999, zero.policy = TRUE)
  z    <- as.numeric(scale(x))
  lagz <- spdep::lag.listw(lw, z, zero.policy = TRUE)
  pval <- lm_res[, ncol(lm_res)]
  padj <- p.adjust(pval, "BH")
  quad <- case_when(z > 0 & lagz > 0 ~ "High-High",
                    z < 0 & lagz < 0 ~ "Low-Low",
                    z > 0 & lagz < 0 ~ "High-Low",
                    TRUE             ~ "Low-High")
  gz <- as.numeric(spdep::localG(x, lw_star, zero.policy = TRUE))
  d %>% select(dcode, district, province) %>%
    mutate(period = p, log_rr = x,
           Ii = lm_res[, "Ii"], p_value = pval, p_fdr = padj,
           cluster     = if_else(pval > 0.05, "Not significant", quad),
           cluster_fdr = if_else(padj > 0.05, "Not significant", quad),
           gi_star = gz,
           hotspot = case_when(2 * pnorm(abs(gz), lower.tail = FALSE) > 0.05 ~ "Not significant",
                               gz > 0 ~ "Hot spot", TRUE ~ "Cold spot"))
})
write.csv(lisa_post, file.path(DIR_RESULTS, "lisa_posterior.csv"), row.names = FALSE)
print(lisa_post %>% count(period, cluster) %>% as.data.frame())

persistent_post <- lisa_post %>%
  filter(cluster == "High-High") %>%
  count(dcode, district, province, name = "n_periods") %>%
  arrange(desc(n_periods), district)
write.csv(persistent_post, file.path(DIR_RESULTS, "lisa_posterior_persistent.csv"),
          row.names = FALSE)
print(as.data.frame(persistent_post))

## Global Moran's I of the posterior surface, period by period
moran_post <- map_dfr(PERIOD_LABELS, function(p) {
  d <- district_rr %>% filter(period == p) %>% arrange(dcode)
  mi <- spdep::moran.mc(log(d$rr), lw, nsim = 9999, zero.policy = TRUE)
  tibble::tibble(period = p, moran_I = unname(mi$statistic), p_value = mi$p.value)
})
write.csv(moran_post, file.path(DIR_RESULTS, "moran_posterior.csv"), row.names = FALSE)
print(as.data.frame(moran_post), digits = 3)

## ===========================================================================
## 4. Model B - individual-level discrete-time survival model
## ===========================================================================

pp <- readRDS(file.path(DIR_DERIVED, "person_period.rds"))

## Birth order and the preceding birth interval overlap exactly on first
## births: birth_order == "1" if and only if bi_cat == "First birth". Entering
## both as dummies makes the pair unidentifiable and inflates the posterior
## variance of each. The interval variable is therefore collapsed to three
## levels among later births, with first births assigned to the reference
## interval; the first-birth effect is then carried by the birth-order term,
## and the interval coefficients are read as the effect of a short interval
## among second and later births.
mb <- pp %>%
  mutate(
    id_space = as.integer(factor(dcode, levels = district_sf$dcode)),
    seg_f    = factor(seg_f, levels = AGE_SEGMENTS$label),
    bi_cat   = factor(if_else(as.character(bi_cat) == "First birth",
                              ">=36 months", as.character(bi_cat)),
                      levels = c(">=36 months", "24-35 months", "<24 months")),
    medu     = relevel(factor(medu),   ref = "Secondary or higher"),
    wealth   = relevel(factor(wealth), ref = "Richest")
  ) %>%
  select(died, seg_f, sex, multiple, birth_order, bi_cat, mage_cat, medu,
         wealth, residence, period, water_imp, sanit_imp, electricity,
         clean_fuel, media, mwork, hh_head_f, ethnicity, alt_km, id_space,
         province, wt, psu, strata, survey_year, age_at_int,
         parity, hhsize) %>%
  tidyr::drop_na(died, seg_f, sex, multiple, birth_order, bi_cat, mage_cat,
                 medu, wealth, residence, period, ethnicity, id_space) %>%
  ## Clusters and strata are numbered afresh in every survey round, so both
  ## identifiers are made unique by prefixing the survey year.
  mutate(id_clust   = as.integer(factor(paste(survey_year, psu))),
         psu_uid    = paste(survey_year, psu),
         strata_uid = paste(survey_year, strata),
         ethnicity  = relevel(factor(ethnicity), ref = "Brahmin/Chhetri"))

## Water, sanitation, electricity and cooking fuel are unrecorded for exactly
## the same children -- those whose mother was not a de jure resident of the
## sampled household -- so their four "Not recorded" dummies are identical
## columns and are not separately identifiable in a regression. Fitted as they
## stand they alias, splitting one effect across four terms with enormous
## posterior variance. They are therefore replaced by a single household-survey
## indicator, leaving the four environmental contrasts to be estimated among
## households that were actually surveyed. Tree ensembles are unaffected by
## this collinearity, so the machine-learning models in 04 keep the three-level
## coding.
collapse_nr <- function(f) {
  f[f == "Not recorded"] <- levels(f)[1]
  droplevels(f)
}
mb <- mb %>%
  mutate(
    hh_env = factor(if_else(as.character(water_imp) == "Not recorded",
                            "Not recorded", "Recorded"),
                    levels = c("Recorded", "Not recorded")),
    across(c(water_imp, sanit_imp, electricity, clean_fuel), collapse_nr)
  )
message("[05] households not surveyed (mother a visitor): ",
        sum(mb$hh_env == "Not recorded"), " person-periods")

message("[05] individual-level model rows: ", nrow(mb), " | deaths: ", sum(mb$died))

## Total children ever born and household size are excluded because the
## outcome itself changes them (see 04); they return only in a sensitivity
## model below. A cluster-level iid effect absorbs the correlation between
## children sampled from the same enumeration area, so the BYM2 district
## variance -- and the median odds ratio derived from it -- measures
## between-district variation net of within-district clustering.
f_ind <- died ~ seg_f + sex + multiple + birth_order + bi_cat + mage_cat +
  medu + wealth + ethnicity + residence + period + water_imp +
  sanit_imp + electricity + clean_fuel + hh_env + media + mwork + hh_head_f +
  scale(alt_km) +
  f(id_space, model = "bym2", graph = g, scale.model = TRUE,
    constr = TRUE, hyper = pc_bym2) +
  f(id_clust, model = "iid", hyper = pc_prec)

message("[05] fitting individual-level model ...")
fit_ind <- inla(f_ind, family = "binomial", data = mb, Ntrials = 1,
                control.predictor = list(compute = FALSE),
                control.compute = list(dic = TRUE, waic = TRUE, config = FALSE),
                control.fixed = list(prec = 0.1, prec.intercept = 0.01),
                control.inla = list(strategy = "adaptive", int.strategy = "eb"),
                verbose = FALSE)

fixed <- as.data.frame(fit_ind$summary.fixed) %>%
  tibble::rownames_to_column("term") %>%
  transmute(term,
            or     = exp(mean),
            or_lo  = exp(`0.025quant`),
            or_hi  = exp(`0.975quant`),
            post_sd = sd,
            prob_gt_1 = 1 - pnorm(0, mean = mean, sd = sd))
write.csv(fixed, file.path(DIR_RESULTS, "inla_fixed_effects.csv"), row.names = FALSE)
print(as.data.frame(fixed %>% filter(!grepl("seg_f", term))), digits = 3)

hyper_ind <- as.data.frame(fit_ind$summary.hyperpar) %>%
  tibble::rownames_to_column("component")
write.csv(hyper_ind, file.path(DIR_RESULTS, "inla_individual_hyper.csv"), row.names = FALSE)
print(hyper_ind, digits = 4)

## Residual district effects after full adjustment ---------------------------
resid_space <- fit_ind$summary.random$id_space[seq_len(nrow(district_sf)), ]
resid_tbl <- tibble::tibble(
  dcode = district_sf$dcode, district = district_sf$district,
  province = district_sf$province,
  resid_log_or = resid_space$mean,
  resid_or = exp(resid_space$mean),
  lo = exp(resid_space$`0.025quant`), hi = exp(resid_space$`0.975quant`)
) %>% arrange(desc(resid_or))
write.csv(resid_tbl, file.path(DIR_RESULTS, "inla_residual_district.csv"), row.names = FALSE)
saveRDS(resid_tbl, file.path(DIR_DERIVED, "inla_residual_district.rds"))
print(head(as.data.frame(resid_tbl), 12), digits = 3)

## Median odds ratio: the median relative increase in the odds of dying that a
## child would experience on moving between two randomly chosen districts.
## Taken over the posterior of each precision rather than its posterior mean,
## so that the median odds ratio carries a credible interval.
mor_from <- function(fit, name) {
  prec <- INLA::inla.rmarginal(4000, fit$marginals.hyperpar[[name]])
  m <- exp(sqrt(2 / prec) * qnorm(0.75))
  c(median = median(m), lo = unname(quantile(m, 0.025)),
    hi = unname(quantile(m, 0.975)))
}
mor_d <- mor_from(fit_ind, "Precision for id_space")
mor_c <- mor_from(fit_ind, "Precision for id_clust")
message("[05] median odds ratio between districts = ", round(mor_d["median"], 3),
        " (", round(mor_d["lo"], 3), "-", round(mor_d["hi"], 3), ")",
        " | between clusters = ", round(mor_c["median"], 3))
write.csv(tibble::tibble(level = c("District", "Cluster"),
                         median_odds_ratio = c(mor_d["median"], mor_c["median"]),
                         lo = c(mor_d["lo"], mor_c["lo"]),
                         hi = c(mor_d["hi"], mor_c["hi"])),
          file.path(DIR_RESULTS, "inla_mor.csv"), row.names = FALSE)

## ===========================================================================
## 5. Design-based sensitivity analysis
## ===========================================================================
## The Bayesian models above are unweighted and include the design variables
## (urban/rural residence, province, wealth) as covariates. As a check, the
## same fixed-effect structure is refitted as a survey-weighted logistic
## regression with Taylor-linearised standard errors that respect the
## stratified cluster design. Cluster and stratum identifiers restart in every
## round, so the survey-prefixed identifiers are used; otherwise cluster 5 of
## 2011 and cluster 5 of 2016 would be treated as one sampling unit.

options(survey.lonely.psu = "adjust")
des <- svydesign(ids = ~psu_uid, strata = ~strata_uid, weights = ~wt,
                 data = mb, nest = TRUE)
f_svy <- died ~ seg_f + sex + multiple + birth_order + bi_cat + mage_cat +
  medu + wealth + ethnicity + residence + period + water_imp +
  sanit_imp + electricity + clean_fuel + hh_env + media + mwork + hh_head_f
svy_fit <- svyglm(f_svy, design = des, family = quasibinomial())
svy_tbl <- broom::tidy(svy_fit, conf.int = TRUE, exponentiate = TRUE)
write.csv(svy_tbl, file.path(DIR_RESULTS, "svy_logistic.csv"), row.names = FALSE)

## ===========================================================================
## 6. Sensitivity analyses for the individual-level model
## ===========================================================================
## (a) Five-year window. Recall of dates and ages at death degrades with time
##     since the event, so the model is refitted on births in the 60 months
##     before interview -- the window DHS uses for national rates.
## (b) Outcome-dependent covariates. Total children ever born and household
##     size are re-entered to show how far the main estimates would move had
##     they been retained.
fit_sens <- function(formula, data, label) {
  message("[05] sensitivity: ", label, " (", nrow(data), " person-periods, ",
          sum(data$died), " deaths)")
  f <- inla(formula, family = "binomial", data = data, Ntrials = 1,
            control.predictor = list(compute = FALSE),
            control.fixed = list(prec = 0.1, prec.intercept = 0.01),
            control.inla = list(strategy = "adaptive", int.strategy = "eb"),
            verbose = FALSE)
  as.data.frame(f$summary.fixed) %>%
    tibble::rownames_to_column("term") %>%
    transmute(analysis = label, term, or = exp(mean),
              or_lo = exp(`0.025quant`), or_hi = exp(`0.975quant`))
}

mb5 <- mb %>%
  filter(age_at_int < 60) %>%
  mutate(id_clust = as.integer(factor(id_clust)), period = droplevels(period))
f_sens_outcome <- update(f_ind, . ~ . + scale(parity) + scale(hhsize))

sens <- bind_rows(
  fixed %>% transmute(analysis = "Main model", term, or, or_lo, or_hi),
  fit_sens(f_ind, mb5, "Births in the 5 years before interview"),
  fit_sens(f_sens_outcome, tidyr::drop_na(mb, parity, hhsize),
           "Adding children ever born and household size")
)
write.csv(sens, file.path(DIR_RESULTS, "inla_sensitivity.csv"), row.names = FALSE)

message("[05] done.")
