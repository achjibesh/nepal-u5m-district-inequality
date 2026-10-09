## ---------------------------------------------------------------------------
## 04_machine_learning.R
##
## Machine-learning analysis of the discrete-time hazard of child death.
##
## Design
##   * Unit of analysis: the person-period (child x DHS age segment). Modelling
##     the segment-specific hazard rather than a binary under-five death
##     indicator removes the right-censoring bias that affects children born
##     less than five years before the interview.
##   * Children (not person-periods) are split 70/30 into a development and a
##     held-out test set, stratified on the child-level outcome and on survey
##     round, so that all person-periods of a child stay in the same partition
##     and no information leaks between them.
##   * Four learners are compared: a penalised logistic regression (elastic
##     net), a random forest, gradient-boosted trees, and a logistic
##     super-learner stacked on the out-of-fold predictions of the other three.
##   * The development set is split into 10 child-grouped folds. Hyperparameter
##     configurations are scored on 5 of them (each fit trained on the other
##     nine tenths); the stack uses out-of-fold predictions from all 10. The
##     test set is touched once.
##   * Spatial transferability is assessed separately by leave-one-province-out
##     cross-validation on the full sample.
##   * The fitted boosting model is interpreted with exact TreeSHAP values,
##     overall and stratified by province.
##
## Outputs (results/): ml_tuning.csv, ml_performance.csv, ml_perf_child.csv,
##   ml_spatial_cv.csv, shap_importance.csv, shap_by_province.csv,
##   pdp_top_features.csv, ml_calibration.csv
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))

suppressPackageStartupMessages({
  library(xgboost)
  library(ranger)
  library(glmnet)
  library(yardstick)
  library(shapviz)
})

pp     <- readRDS(file.path(DIR_DERIVED, "person_period.rds"))
births <- readRDS(file.path(DIR_DERIVED, "births.rds"))

## ===========================================================================
## 1. Feature set
## ===========================================================================
## Predictors are grouped into four blocks that map onto the Mosley-Chen
## framework for child survival: the age segment itself (the baseline hazard),
## child and birth characteristics, maternal characteristics, and household /
## environmental / geographic context.
##
## Two rules govern the list.
##  * One encoding per concept. Entering both a grouped and a continuous
##    version of the same quantity (birth order, birth interval, maternal age,
##    education, wealth) splits its SHAP importance across two columns and
##    understates it. The continuous form is kept; tree ensembles find their
##    own cut-points.
##  * Nothing the outcome itself can change. Total children ever born (v201)
##    rises when a dead child is replaced by a later birth, and household size
##    (v136) falls when a child dies because the dead are not household
##    members. Both are recorded at interview, after the death, so each partly
##    encodes the outcome. Neither is used; birth order, fixed at the child's
##    birth, carries the family-size signal instead.

FEATURES <- list(
  hazard    = c("seg_f"),
  child     = c("sex", "multiple", "bord_num", "bi_num"),
  maternal  = c("mage_birth", "medu_years", "mwork", "media", "ethnicity"),
  household = c("wealth_score", "residence", "hh_head_f",
                "water_imp", "sanit_imp", "electricity", "clean_fuel"),
  context   = c("province", "alt_km", "period", "round")
)
FEAT <- unlist(FEATURES, use.names = FALSE)

FEATURE_LABELS <- c(
  seg_f = "Age segment", sex = "Sex of child", multiple = "Multiple birth",
  bord_num = "Birth order", bi_num = "Preceding birth interval (months)",
  mage_birth = "Maternal age at birth (years)",
  medu_years = "Maternal years of schooling", mwork = "Maternal employment",
  media = "Maternal media exposure", ethnicity = "Caste/ethnicity",
  wealth_score = "Household wealth index score", residence = "Place of residence",
  hh_head_f = "Sex of household head",
  water_imp = "Improved drinking water", sanit_imp = "Improved sanitation",
  electricity = "Household electricity", clean_fuel = "Clean cooking fuel",
  province = "Province", alt_km = "Altitude (km)",
  period = "Birth period", round = "Survey round"
)

dat <- pp %>%
  select(child_id, died, expo, all_of(FEAT), dcode, district, survey_year,
         u5_death, complete_exposure, wt, psu) %>%
  mutate(
    ## a missing preceding birth interval means a first birth; it is coded
    ## explicitly rather than imputed
    bi_num = tidyr::replace_na(bi_num, -1),
    across(where(is.character), as.factor)
  )

## Remaining missingness is negligible; complete cases are used.
na_by_feat <- sapply(dat[FEAT], function(x) mean(is.na(x)))
message("[04] features with any missingness: ",
        paste(names(na_by_feat)[na_by_feat > 0], collapse = ", "))
dat <- dat %>% filter(if_all(all_of(FEAT), ~ !is.na(.)))
message("[04] person-periods: ", nrow(dat), " | segment deaths: ", sum(dat$died))

## One-hot design matrix, shared by every learner ---------------------------
X_all <- model.matrix(~ . - 1, data = dat[, FEAT, drop = FALSE])
colnames(X_all) <- make.names(colnames(X_all))
y_all <- dat$died
message("[04] design matrix: ", nrow(X_all), " x ", ncol(X_all))

## ===========================================================================
## 2. Development / test partition (grouped by child)
## ===========================================================================

child_key <- dat %>%
  group_by(child_id) %>%
  summarise(u5 = max(died), round = first(survey_year), .groups = "drop")

set.seed(SEED)
child_key <- child_key %>%
  group_by(u5, round) %>%
  mutate(part = if_else(runif(n()) < 0.70, "dev", "test")) %>%
  ungroup()

dev_ids  <- child_key$child_id[child_key$part == "dev"]
test_ids <- child_key$child_id[child_key$part == "test"]
idx_dev  <- which(dat$child_id %in% dev_ids)
idx_test <- which(dat$child_id %in% test_ids)

message("[04] development children: ", length(dev_ids),
        " (deaths ", sum(child_key$u5[child_key$part == "dev"]), ")",
        " | test children: ", length(test_ids),
        " (deaths ", sum(child_key$u5[child_key$part == "test"]), ")")

## 10 grouped folds inside the development set ------------------------------
K <- 10
set.seed(SEED + 1)
fold_of_child <- setNames(sample(rep_len(1:K, length(dev_ids))), dev_ids)
fold_dev <- unname(fold_of_child[as.character(dat$child_id[idx_dev])])

## ===========================================================================
## 3. Learners
## ===========================================================================

fit_xgb <- function(Xtr, ytr, par, nrounds) {
  xgboost::xgboost(
    data = Xtr, label = ytr, nrounds = nrounds,
    max_depth = par$max_depth, eta = par$eta,
    subsample = par$subsample, colsample_bytree = par$colsample,
    min_child_weight = par$min_child_weight,
    scale_pos_weight = par$scale_pos_weight,
    objective = "binary:logistic", eval_metric = "logloss",
    tree_method = "hist", nthread = 4, verbose = 0)
}

fit_rf <- function(Xtr, ytr, par) {
  ranger::ranger(x = Xtr, y = factor(ytr, levels = c(0, 1)), probability = TRUE,
                 num.trees = par$num_trees, mtry = par$mtry,
                 min.node.size = par$min_node, max.depth = par$max_depth,
                 sample.fraction = par$sample_fraction,
                 class.weights = c(1, par$pos_weight),
                 num.threads = 4, seed = SEED, verbose = FALSE)
}

fit_enet <- function(Xtr, ytr, alpha = 0.5) {
  cvfit <- glmnet::cv.glmnet(Xtr, ytr, family = "binomial", alpha = alpha,
                             nfolds = 5, type.measure = "deviance")
  list(fit = cvfit, lambda = cvfit$lambda.min)
}

pred_of <- function(model, kind, Xnew) {
  switch(kind,
         xgb  = predict(model, Xnew),
         rf   = predict(model, Xnew, num.threads = 4)$predictions[, "1"],
         enet = as.numeric(predict(model$fit, Xnew, s = model$lambda,
                                   type = "response")))
}

## ===========================================================================
## 4. Hyperparameter search (grouped CV inside the development set: each
##    configuration is scored on folds 1-5 of the 10, to halve the run time)
## ===========================================================================

pr_auc_vec <- function(truth, prob) {
  d <- tibble::tibble(truth = factor(truth, levels = c(1, 0)), prob = prob)
  yardstick::pr_auc_vec(d$truth, d$prob, event_level = "first")
}
roc_auc_vec2 <- function(truth, prob) {
  d <- tibble::tibble(truth = factor(truth, levels = c(1, 0)), prob = prob)
  yardstick::roc_auc_vec(d$truth, d$prob, event_level = "first")
}

cv_score <- function(kind, par, folds_used = 1:K) {
  oof <- rep(NA_real_, length(idx_dev))
  for (k in folds_used) {
    tr <- idx_dev[fold_dev != k]; va <- idx_dev[fold_dev == k]
    m  <- switch(kind,
                 xgb  = fit_xgb(X_all[tr, ], y_all[tr], par, par$nrounds),
                 rf   = fit_rf(X_all[tr, ], y_all[tr], par),
                 enet = fit_enet(X_all[tr, ], y_all[tr], par$alpha))
    oof[fold_dev == k] <- pred_of(m, kind, X_all[va, ])
  }
  keep <- !is.na(oof)
  list(oof = oof,
       roc = roc_auc_vec2(y_all[idx_dev][keep], oof[keep]),
       pr  = pr_auc_vec(y_all[idx_dev][keep], oof[keep]))
}

set.seed(SEED + 2)
N_CONFIG <- 12
grid_xgb <- tibble::tibble(
  max_depth        = sample(2:6, N_CONFIG, TRUE),
  eta              = round(exp(runif(N_CONFIG, log(0.01), log(0.20))), 3),
  subsample        = round(runif(N_CONFIG, 0.6, 1.0), 2),
  colsample        = round(runif(N_CONFIG, 0.5, 1.0), 2),
  min_child_weight = sample(c(1, 5, 10, 25), N_CONFIG, TRUE),
  scale_pos_weight = sample(c(1, 5, 20), N_CONFIG, TRUE),
  nrounds          = sample(c(200, 400, 600), N_CONFIG, TRUE)
)
grid_rf <- tibble::tibble(
  num_trees       = 500,
  mtry            = sample(3:15, 6, TRUE),
  min_node        = sample(c(10, 25, 50, 100), 6, TRUE),
  max_depth       = sample(c(0, 6, 10, 15), 6, TRUE),
  sample_fraction = round(runif(6, 0.6, 1.0), 2),
  pos_weight      = sample(c(1, 5, 20), 6, TRUE)
)

## The search is the expensive part of the pipeline (roughly 45 minutes on four
## threads). Its result is cached: if ml_tuning.csv already exists the selected
## configurations are read back, unless U5M_RETUNE=1 forces a fresh search.
TUNE_FILE <- file.path(DIR_RESULTS, "ml_tuning.csv")
if (file.exists(TUNE_FILE) && Sys.getenv("U5M_RETUNE") != "1") {
  message("[04] reusing cached hyperparameter search (", TUNE_FILE, ")")
  tuned <- read.csv(TUNE_FILE)
  num   <- function(x) suppressWarnings(as.numeric(as.character(x)))
  tx <- tuned[tuned$learner == "Gradient boosting", ]
  tr <- tuned[tuned$learner == "Random forest", ]
  best_xgb <- lapply(tx[which.max(num(tx$pr_auc)), names(grid_xgb)], num)
  best_rf  <- lapply(tr[which.max(num(tr$pr_auc)), names(grid_rf)],  num)
} else {
  message("[04] tuning gradient boosting ...")
  tune_xgb <- map_dfr(seq_len(nrow(grid_xgb)), function(i) {
    s <- cv_score("xgb", as.list(grid_xgb[i, ]), folds_used = 1:5)
    bind_cols(grid_xgb[i, ], tibble::tibble(learner = "Gradient boosting",
                                            roc_auc = s$roc, pr_auc = s$pr))
  })
  message("[04] tuning random forest ...")
  tune_rf <- map_dfr(seq_len(nrow(grid_rf)), function(i) {
    s <- cv_score("rf", as.list(grid_rf[i, ]), folds_used = 1:5)
    bind_cols(grid_rf[i, ], tibble::tibble(learner = "Random forest",
                                           roc_auc = s$roc, pr_auc = s$pr))
  })
  best_xgb <- as.list(tune_xgb[which.max(tune_xgb$pr_auc), names(grid_xgb)])
  best_rf  <- as.list(tune_rf[which.max(tune_rf$pr_auc),  names(grid_rf)])
  write.csv(bind_rows(tune_xgb %>% mutate(across(everything(), as.character)),
                      tune_rf  %>% mutate(across(everything(), as.character))),
            TUNE_FILE, row.names = FALSE)
}
message("[04] best boosting config: ", paste(names(best_xgb), unlist(best_xgb),
                                             sep = "=", collapse = ", "))
message("[04] best forest config: ", paste(names(best_rf), unlist(best_rf),
                                           sep = "=", collapse = ", "))

## ===========================================================================
## 5. Out-of-fold predictions, stacking, and held-out evaluation
## ===========================================================================

message("[04] computing out-of-fold predictions for the stack ...")
oof_xgb  <- cv_score("xgb",  best_xgb,               1:K)$oof
oof_rf   <- cv_score("rf",   best_rf,                1:K)$oof
oof_enet <- cv_score("enet", list(alpha = 0.5),      1:K)$oof

stack_df <- data.frame(y = y_all[idx_dev],
                       xgb = qlogis(pmin(pmax(oof_xgb, 1e-6), 1 - 1e-6)),
                       rf  = qlogis(pmin(pmax(oof_rf,  1e-6), 1 - 1e-6)),
                       enet = qlogis(pmin(pmax(oof_enet, 1e-6), 1 - 1e-6)))
stack_fit <- glm(y ~ xgb + rf + enet, data = stack_df, family = binomial())

## refit each learner on the whole development set --------------------------
final_xgb  <- fit_xgb(X_all[idx_dev, ], y_all[idx_dev], best_xgb, best_xgb$nrounds)
final_rf   <- fit_rf(X_all[idx_dev, ],  y_all[idx_dev], best_rf)
final_enet <- fit_enet(X_all[idx_dev, ], y_all[idx_dev], 0.5)

test_pred <- tibble::tibble(
  child_id = dat$child_id[idx_test],
  seg      = dat$seg_f[idx_test],
  y        = y_all[idx_test],
  `Gradient boosting` = pred_of(final_xgb,  "xgb",  X_all[idx_test, ]),
  `Random forest`     = pred_of(final_rf,   "rf",   X_all[idx_test, ]),
  `Elastic-net logistic regression` = pred_of(final_enet, "enet", X_all[idx_test, ])
)
test_pred$`Stacked ensemble` <- predict(
  stack_fit, newdata = data.frame(
    xgb  = qlogis(pmin(pmax(test_pred$`Gradient boosting`, 1e-6), 1 - 1e-6)),
    rf   = qlogis(pmin(pmax(test_pred$`Random forest`, 1e-6), 1 - 1e-6)),
    enet = qlogis(pmin(pmax(test_pred$`Elastic-net logistic regression`, 1e-6), 1 - 1e-6))),
  type = "response")

LEARNERS <- c("Elastic-net logistic regression", "Random forest",
              "Gradient boosting", "Stacked ensemble")

## Bootstrap confidence intervals, resampling children ----------------------
boot_metric <- function(truth, prob, groups, fun, R = 500, seed = SEED) {
  set.seed(seed)
  g   <- unique(groups)
  idx <- split(seq_along(groups), groups)
  out <- numeric(R)
  for (r in seq_len(R)) {
    take <- unlist(idx[sample(g, length(g), replace = TRUE)], use.names = FALSE)
    out[r] <- tryCatch(fun(truth[take], prob[take]), error = function(e) NA_real_)
  }
  stats::quantile(out, c(0.025, 0.975), na.rm = TRUE)
}

brier <- function(truth, prob) mean((prob - truth)^2)

perf <- map_dfr(LEARNERS, function(L) {
  p <- test_pred[[L]]; y <- test_pred$y; g <- test_pred$child_id
  ci_roc <- boot_metric(y, p, g, roc_auc_vec2, R = 300)
  ci_pr  <- boot_metric(y, p, g, pr_auc_vec,   R = 300)
  tibble::tibble(
    learner = L, level = "Person-period (segment hazard)",
    roc_auc = roc_auc_vec2(y, p), roc_lo = ci_roc[1], roc_hi = ci_roc[2],
    pr_auc  = pr_auc_vec(y, p),   pr_lo  = ci_pr[1],  pr_hi  = ci_pr[2],
    brier   = brier(y, p),
    prevalence = mean(y))
})

## Written out as soon as they exist, so that a failure further down the script
## does not discard an hour of model fitting.
write.csv(perf,    file.path(DIR_RESULTS, "ml_performance.csv"), row.names = FALSE)
saveRDS(test_pred, file.path(DIR_DERIVED, "test_predictions.rds"))
saveRDS(list(best_xgb = best_xgb, best_rf = best_rf,
             stack = coef(stack_fit)),
        file.path(DIR_DERIVED, "ml_config.rds"))

## ---- child-level evaluation ----------------------------------------------
## The person-period hazards of a child are combined into a predicted
## probability of dying before age five, 1 - prod(1 - h_s), evaluated over all
## eight age segments. Discrimination is then assessed among children whose
## survival to age five was actually observed.

full_grid <- dat[idx_test, ] %>%
  distinct(child_id, .keep_all = TRUE) %>%
  select(-seg_f) %>%
  tidyr::crossing(seg_f = factor(AGE_SEGMENTS$label, levels = AGE_SEGMENTS$label))

## The grid can lack factor levels that occur in the full sample, so the design
## matrix is aligned to the training columns rather than subset directly.
align_cols <- function(X, template) {
  colnames(X) <- make.names(colnames(X))
  out <- matrix(0, nrow(X), length(template),
                dimnames = list(NULL, template))
  shared <- intersect(colnames(X), template)
  out[, shared] <- X[, shared, drop = FALSE]
  out
}
Xg <- align_cols(model.matrix(~ . - 1, data = full_grid[, FEAT, drop = FALSE]),
                 colnames(X_all))

child_pred <- map_dfr(LEARNERS[LEARNERS != "Stacked ensemble"], function(L) {
  kind <- switch(L, `Gradient boosting` = "xgb", `Random forest` = "rf",
                 `Elastic-net logistic regression` = "enet")
  mod  <- switch(kind, xgb = final_xgb, rf = final_rf, enet = final_enet)
  h    <- pred_of(mod, kind, Xg)
  full_grid %>%
    mutate(h = h) %>%
    group_by(child_id) %>%
    summarise(q5 = 1 - prod(1 - h),
              u5_death = first(u5_death),
              complete_exposure = first(complete_exposure), .groups = "drop") %>%
    mutate(learner = L)
})

perf_child <- child_pred %>%
  filter(complete_exposure) %>%
  group_by(learner) %>%
  group_modify(~ tibble::tibble(
    level      = "Child (probability of dying before age five)",
    roc_auc    = roc_auc_vec2(.x$u5_death, .x$q5),
    pr_auc     = pr_auc_vec(.x$u5_death, .x$q5),
    brier      = brier(.x$u5_death, .x$q5),
    prevalence = mean(.x$u5_death),
    n          = nrow(.x))) %>%
  ungroup()

print(as.data.frame(perf %>% select(learner, roc_auc, roc_lo, roc_hi, pr_auc, brier)))
print(as.data.frame(perf_child))
write.csv(perf,       file.path(DIR_RESULTS, "ml_performance.csv"), row.names = FALSE)
write.csv(perf_child, file.path(DIR_RESULTS, "ml_perf_child.csv"),  row.names = FALSE)
saveRDS(test_pred,    file.path(DIR_DERIVED, "test_predictions.rds"))
saveRDS(child_pred,   file.path(DIR_DERIVED, "child_predictions.rds"))

## ---- calibration ----------------------------------------------------------
calib <- map_dfr(LEARNERS, function(L) {
  p <- test_pred[[L]]
  tibble::tibble(learner = L, p = p, y = test_pred$y) %>%
    mutate(bin = ntile(p, 10)) %>%
    group_by(learner, bin) %>%
    summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop")
})
calib_slope <- map_dfr(LEARNERS, function(L) {
  p  <- pmin(pmax(test_pred[[L]], 1e-6), 1 - 1e-6)
  m  <- glm(test_pred$y ~ qlogis(p), family = binomial())
  tibble::tibble(learner = L, calib_intercept = coef(m)[1], calib_slope = coef(m)[2])
})
write.csv(calib,       file.path(DIR_RESULTS, "ml_calibration.csv"), row.names = FALSE)
write.csv(calib_slope, file.path(DIR_RESULTS, "ml_calibration_slope.csv"), row.names = FALSE)

## ===========================================================================
## 6. Spatial transferability: leave-one-province-out cross-validation
## ===========================================================================
## Each province is held out in turn and the model is trained on the other
## six. This asks whether a risk model learned in one part of Nepal transfers
## to a region it has never seen, which random cross-validation cannot answer
## because neighbouring clusters leak into both partitions.

provs <- levels(droplevels(dat$province))
spatial_cv <- map_dfr(provs, function(pv) {
  tr <- which(dat$province != pv); te <- which(dat$province == pv)
  if (sum(y_all[te]) < 5) return(NULL)
  m_x <- fit_xgb(X_all[tr, ], y_all[tr], best_xgb, best_xgb$nrounds)
  m_e <- fit_enet(X_all[tr, ], y_all[tr], 0.5)
  tibble::tibble(
    province = pv, n_person_periods = length(te), deaths = sum(y_all[te]),
    roc_xgb  = roc_auc_vec2(y_all[te], pred_of(m_x, "xgb",  X_all[te, ])),
    pr_xgb   = pr_auc_vec(  y_all[te], pred_of(m_x, "xgb",  X_all[te, ])),
    roc_enet = roc_auc_vec2(y_all[te], pred_of(m_e, "enet", X_all[te, ])),
    pr_enet  = pr_auc_vec(  y_all[te], pred_of(m_e, "enet", X_all[te, ])))
})
print(as.data.frame(spatial_cv))
write.csv(spatial_cv, file.path(DIR_RESULTS, "ml_spatial_cv.csv"), row.names = FALSE)

## ===========================================================================
## 7. Interpretation: TreeSHAP
## ===========================================================================
## The boosting model is refitted on the complete sample and explained with
## exact TreeSHAP contributions on the log-odds scale. Contributions from the
## dummy columns of one categorical variable are summed back to the variable.

final_full <- fit_xgb(X_all, y_all, best_xgb, best_xgb$nrounds)
saveRDS(list(model = xgboost::xgb.save.raw(final_full), best = best_xgb,
             cols = colnames(X_all)),
        file.path(DIR_DERIVED, "xgb_final.rds"))

set.seed(SEED)
shap_idx <- sample(seq_len(nrow(X_all)), min(60000, nrow(X_all)))
sv <- shapviz::shapviz(final_full, X_pred = X_all[shap_idx, ], X = X_all[shap_idx, ])
S  <- shapviz::get_shap_values(sv)

## map one-hot columns back to their source variable
col_to_var <- function(cols, feats) {
  vapply(cols, function(cn) {
    hit <- feats[startsWith(cn, make.names(feats))]
    if (length(hit) == 0) return(NA_character_)
    hit[which.max(nchar(hit))]
  }, character(1))
}
map_var <- col_to_var(colnames(S), FEAT)

shap_importance <- tibble::tibble(column = colnames(S),
                                  variable = map_var,
                                  mean_abs = colMeans(abs(S))) %>%
  group_by(variable) %>%
  summarise(mean_abs_shap = sum(mean_abs), .groups = "drop") %>%
  arrange(desc(mean_abs_shap)) %>%
  mutate(label = dplyr::coalesce(FEATURE_LABELS[variable], variable),
         block = case_when(variable %in% FEATURES$hazard    ~ "Baseline hazard",
                           variable %in% FEATURES$child     ~ "Child and birth",
                           variable %in% FEATURES$maternal  ~ "Maternal",
                           variable %in% FEATURES$household ~ "Household",
                           TRUE                             ~ "Geographic context"),
         rel = mean_abs_shap / sum(mean_abs_shap))
print(as.data.frame(shap_importance), digits = 3)
write.csv(shap_importance, file.path(DIR_RESULTS, "shap_importance.csv"), row.names = FALSE)

## Province-stratified SHAP: does the same risk factor carry the same weight
## everywhere in Nepal?
prov_of_row <- dat$province[shap_idx]
shap_by_province <- map_dfr(levels(droplevels(prov_of_row)), function(pv) {
  rows <- which(prov_of_row == pv)
  tibble::tibble(province = pv, column = colnames(S),
                 variable = map_var,
                 mean_abs = colMeans(abs(S[rows, , drop = FALSE]))) %>%
    group_by(province, variable) %>%
    summarise(mean_abs_shap = sum(mean_abs), .groups = "drop") %>%
    mutate(rel = mean_abs_shap / sum(mean_abs_shap))
})
write.csv(shap_by_province, file.path(DIR_RESULTS, "shap_by_province.csv"), row.names = FALSE)

## Long SHAP table for the beeswarm figure ----------------------------------
top_vars <- shap_importance$variable[1:12]
shap_long <- as.data.frame(S) %>%
  mutate(.row = row_number()) %>%
  tidyr::pivot_longer(-.row, names_to = "column", values_to = "shap") %>%
  mutate(variable = map_var[column]) %>%
  filter(variable %in% top_vars) %>%
  left_join(as.data.frame(X_all[shap_idx, ]) %>%
              mutate(.row = row_number()) %>%
              tidyr::pivot_longer(-.row, names_to = "column", values_to = "x"),
            by = c(".row", "column"))
saveRDS(shap_long, file.path(DIR_DERIVED, "shap_long.rds"))

## Partial dependence for the leading continuous predictors -----------------
pdp_one <- function(col, grid_vals) {
  map_dfr(grid_vals, function(g) {
    Xtmp <- X_all[shap_idx, , drop = FALSE]
    Xtmp[, col] <- g
    tibble::tibble(column = col, value = g,
                   yhat = mean(predict(final_full, Xtmp)))
  })
}
cont_cols <- intersect(c("bi_num", "mage_birth", "medu_years", "wealth_score",
                         "alt_km", "bord_num"),
                       colnames(X_all))
pdp <- map_dfr(cont_cols, function(cl) {
  v <- X_all[, cl]
  gr <- unique(round(stats::quantile(v[v > -1], probs = seq(0.02, 0.98, length.out = 20),
                                     na.rm = TRUE), 3))
  pdp_one(cl, gr)
}) %>%
  mutate(variable = column,
         label = dplyr::coalesce(FEATURE_LABELS[variable], variable))
write.csv(pdp, file.path(DIR_RESULTS, "pdp_top_features.csv"), row.names = FALSE)

message("[04] done.")
