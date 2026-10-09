## ---------------------------------------------------------------------------
## 04b_recalibration.R
##
## The hyperparameter search in 04 selects the learners on the area under the
## precision-recall curve, and for both tree ensembles it chooses a large
## positive-class weight. Class weighting improves the ranking of rare events
## but multiplies the predicted odds by the weight, so the raw probabilities
## from those models are badly calibrated and their Brier scores and implied
## absolute risks are not interpretable.
##
## This script therefore recalibrates. Platt scaling (a logistic regression of
## the outcome on the predicted log-odds) is fitted on *development-set
## out-of-fold* predictions and applied to the held-out test set, so no test
## information is used to fit the correction. Recalibration is a monotone
## transformation and so leaves AUROC and AUPRC unchanged; it repairs the Brier
## score, the calibration curve and the child-level probability of dying before
## age five, all of which are reported in the manuscript.
##
## Overwrites: ml_performance.csv, ml_perf_child.csv, ml_calibration.csv,
##             ml_calibration_slope.csv, data_derived/test_predictions.rds
## Writes:     ml_recalibration.csv (the fitted scaling coefficients)
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(xgboost); library(ranger); library(glmnet)
  library(yardstick)
})

pp <- readRDS(file.path(DIR_DERIVED, "person_period.rds"))

## ---- rebuild exactly the design matrix and split used in 04 ---------------
## must match 04_machine_learning.R exactly
FEATURES <- list(
  hazard    = c("seg_f"),
  child     = c("sex", "multiple", "bord_num", "bi_num"),
  maternal  = c("mage_birth", "medu_years", "mwork", "media", "ethnicity"),
  household = c("wealth_score", "residence", "hh_head_f",
                "water_imp", "sanit_imp", "electricity", "clean_fuel"),
  context   = c("province", "alt_km", "period", "round")
)
FEAT <- unlist(FEATURES, use.names = FALSE)

dat <- pp %>%
  select(child_id, died, expo, all_of(FEAT), dcode, district, survey_year,
         u5_death, complete_exposure, wt, psu) %>%
  mutate(bi_num = tidyr::replace_na(bi_num, -1),
         across(where(is.character), as.factor)) %>%
  filter(if_all(all_of(FEAT), ~ !is.na(.)))

X_all <- model.matrix(~ . - 1, data = dat[, FEAT, drop = FALSE])
colnames(X_all) <- make.names(colnames(X_all))
y_all <- dat$died

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

K <- 10
set.seed(SEED + 1)
fold_of_child <- setNames(sample(rep_len(1:K, length(dev_ids))), dev_ids)
fold_dev <- unname(fold_of_child[as.character(dat$child_id[idx_dev])])

stopifnot(identical(length(idx_test),
                    nrow(readRDS(file.path(DIR_DERIVED, "test_predictions.rds")))))

## ---- best hyperparameters, as chosen in 04 --------------------------------
tune <- read.csv(file.path(DIR_RESULTS, "ml_tuning.csv"))
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
tx <- tune[tune$learner == "Gradient boosting", ]
tr <- tune[tune$learner == "Random forest", ]
best_xgb <- as.list(lapply(tx[which.max(num(tx$pr_auc)),
                              c("max_depth","eta","subsample","colsample",
                                "min_child_weight","scale_pos_weight","nrounds")], num))
best_rf <- as.list(lapply(tr[which.max(num(tr$pr_auc)),
                             c("num_trees","mtry","min_node","max_depth",
                               "sample_fraction","pos_weight")], num))

fit_xgb <- function(Xtr, ytr, par) {
  xgboost::xgboost(data = Xtr, label = ytr, nrounds = par$nrounds,
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
lo <- function(p) qlogis(pmin(pmax(p, 1e-7), 1 - 1e-7))

## ---- development-set out-of-fold predictions for the scaling fit ----------
## Five folds are enough to estimate a two-parameter correction.
message("[04b] out-of-fold predictions on the development set ...")
oof <- data.frame(y = y_all[idx_dev], xgb = NA_real_, rf = NA_real_,
                  enet = NA_real_)
for (k in 1:5) {
  tr_i <- idx_dev[fold_dev != k & fold_dev <= 5]
  va_i <- idx_dev[fold_dev == k]
  mx <- fit_xgb(X_all[tr_i, ], y_all[tr_i], best_xgb)
  mr <- fit_rf(X_all[tr_i, ], y_all[tr_i], best_rf)
  me <- glmnet::cv.glmnet(X_all[tr_i, ], y_all[tr_i], family = "binomial",
                          alpha = 0.5, nfolds = 5)
  oof$xgb[fold_dev == k]  <- predict(mx, X_all[va_i, ])
  oof$rf[fold_dev == k]   <- predict(mr, X_all[va_i, ], num.threads = 4)$predictions[, "1"]
  oof$enet[fold_dev == k] <- as.numeric(predict(me, X_all[va_i, ],
                                                s = me$lambda.min, type = "response"))
  message("   fold ", k, " done")
}
oof <- oof[!is.na(oof$xgb), ]

platt_xgb <- glm(y ~ lo(xgb), data = oof, family = binomial())
platt_rf  <- glm(y ~ lo(rf),  data = oof, family = binomial())
recal <- tibble::tibble(
  learner = c("Gradient boosting", "Random forest"),
  intercept = c(coef(platt_xgb)[1], coef(platt_rf)[1]),
  slope     = c(coef(platt_xgb)[2], coef(platt_rf)[2]),
  weight_used = c(best_xgb$scale_pos_weight, best_rf$pos_weight))
print(as.data.frame(recal))
write.csv(recal, file.path(DIR_RESULTS, "ml_recalibration.csv"), row.names = FALSE)

scale_p <- function(p, mod) plogis(coef(mod)[1] + coef(mod)[2] * lo(p))

## ---- apply to the held-out test set ---------------------------------------
test_pred <- readRDS(file.path(DIR_DERIVED, "test_predictions.rds"))
test_pred$`Gradient boosting` <- scale_p(test_pred$`Gradient boosting`, platt_xgb)
test_pred$`Random forest`     <- scale_p(test_pred$`Random forest`, platt_rf)

## The super-learner is refitted on the *recalibrated* development-set
## out-of-fold predictions, so that the ensemble is a logistic stack of the
## same three calibrated learners that the table reports individually.
stack_dev <- data.frame(
  y    = oof$y,
  xgb  = lo(scale_p(oof$xgb, platt_xgb)),
  rf   = lo(scale_p(oof$rf,  platt_rf)),
  enet = lo(oof$enet))
stack_fit <- glm(y ~ xgb + rf + enet, data = stack_dev, family = binomial())
test_pred$`Stacked ensemble` <- as.numeric(predict(
  stack_fit,
  newdata = data.frame(xgb  = lo(test_pred$`Gradient boosting`),
                       rf   = lo(test_pred$`Random forest`),
                       enet = lo(test_pred$`Elastic-net logistic regression`)),
  type = "response"))
write.csv(broom::tidy(stack_fit), file.path(DIR_RESULTS, "ml_stack_weights.csv"),
          row.names = FALSE)

LEARNERS <- c("Elastic-net logistic regression", "Random forest",
              "Gradient boosting", "Stacked ensemble")
pr_auc_v  <- function(t, p) yardstick::pr_auc_vec(factor(t, levels = c(1, 0)), p,
                                                  event_level = "first")
roc_auc_v <- function(t, p) yardstick::roc_auc_vec(factor(t, levels = c(1, 0)), p,
                                                   event_level = "first")
brier <- function(t, p) mean((p - t)^2)

boot_metric <- function(truth, prob, groups, fun, R = 300, seed = SEED) {
  set.seed(seed); g <- unique(groups); idx <- split(seq_along(groups), groups)
  out <- numeric(R)
  for (r in seq_len(R)) {
    take <- unlist(idx[sample(g, length(g), replace = TRUE)], use.names = FALSE)
    out[r] <- tryCatch(fun(truth[take], prob[take]), error = function(e) NA_real_)
  }
  stats::quantile(out, c(0.025, 0.975), na.rm = TRUE)
}

perf <- purrr::map_dfr(LEARNERS, function(L) {
  p <- test_pred[[L]]; y <- test_pred$y; g <- test_pred$child_id
  cr <- boot_metric(y, p, g, roc_auc_v); cp <- boot_metric(y, p, g, pr_auc_v)
  tibble::tibble(learner = L, level = "Person-period (segment hazard)",
                 roc_auc = roc_auc_v(y, p), roc_lo = cr[1], roc_hi = cr[2],
                 pr_auc = pr_auc_v(y, p), pr_lo = cp[1], pr_hi = cp[2],
                 brier = brier(y, p), prevalence = mean(y))
})
print(as.data.frame(perf %>% select(learner, roc_auc, pr_auc, brier)))
write.csv(perf, file.path(DIR_RESULTS, "ml_performance.csv"), row.names = FALSE)
saveRDS(test_pred, file.path(DIR_DERIVED, "test_predictions.rds"))

calib <- purrr::map_dfr(LEARNERS, function(L) {
  tibble::tibble(learner = L, p = test_pred[[L]], y = test_pred$y) %>%
    mutate(bin = ntile(p, 10)) %>%
    group_by(learner, bin) %>%
    summarise(pred = mean(p), obs = mean(y), n = n(), .groups = "drop")
})
calib_slope <- purrr::map_dfr(LEARNERS, function(L) {
  m <- glm(test_pred$y ~ lo(test_pred[[L]]), family = binomial())
  tibble::tibble(learner = L, calib_intercept = coef(m)[1], calib_slope = coef(m)[2])
})
write.csv(calib, file.path(DIR_RESULTS, "ml_calibration.csv"), row.names = FALSE)
write.csv(calib_slope, file.path(DIR_RESULTS, "ml_calibration_slope.csv"), row.names = FALSE)
print(as.data.frame(calib_slope))

## ---- child-level probability of dying before age five ---------------------
message("[04b] recomputing child-level probabilities ...")
final_xgb  <- fit_xgb(X_all[idx_dev, ], y_all[idx_dev], best_xgb)
final_rf   <- fit_rf(X_all[idx_dev, ],  y_all[idx_dev], best_rf)
cvfit      <- glmnet::cv.glmnet(X_all[idx_dev, ], y_all[idx_dev],
                                family = "binomial", alpha = 0.5, nfolds = 5)

grid <- dat[idx_test, ] %>%
  distinct(child_id, .keep_all = TRUE) %>%
  select(-seg_f) %>%
  tidyr::crossing(seg_f = factor(AGE_SEGMENTS$label, levels = AGE_SEGMENTS$label))
align_cols <- function(X, template) {
  colnames(X) <- make.names(colnames(X))
  out <- matrix(0, nrow(X), length(template), dimnames = list(NULL, template))
  shared <- intersect(colnames(X), template)
  out[, shared] <- X[, shared, drop = FALSE]
  out
}
Xg <- align_cols(model.matrix(~ . - 1, data = grid[, FEAT, drop = FALSE]),
                 colnames(X_all))

haz <- list(
  `Gradient boosting` = scale_p(predict(final_xgb, Xg), platt_xgb),
  `Random forest`     = scale_p(predict(final_rf, Xg, num.threads = 4)$predictions[, "1"],
                                platt_rf),
  `Elastic-net logistic regression` =
    as.numeric(predict(cvfit, Xg, s = cvfit$lambda.min, type = "response"))
)

child_pred <- purrr::map_dfr(names(haz), function(L) {
  grid %>% mutate(h = haz[[L]]) %>%
    group_by(child_id) %>%
    summarise(q5 = 1 - prod(1 - h), u5_death = first(u5_death),
              complete_exposure = first(complete_exposure), .groups = "drop") %>%
    mutate(learner = L)
})
perf_child <- child_pred %>%
  filter(complete_exposure) %>%
  group_by(learner) %>%
  group_modify(~ tibble::tibble(
    level = "Child (probability of dying before age five)",
    roc_auc = roc_auc_v(.x$u5_death, .x$q5),
    pr_auc  = pr_auc_v(.x$u5_death, .x$q5),
    brier   = brier(.x$u5_death, .x$q5),
    mean_predicted = mean(.x$q5), observed = mean(.x$u5_death),
    prevalence = mean(.x$u5_death), n = nrow(.x))) %>%
  ungroup()
print(as.data.frame(perf_child))
write.csv(perf_child, file.path(DIR_RESULTS, "ml_perf_child.csv"), row.names = FALSE)
saveRDS(child_pred, file.path(DIR_DERIVED, "child_predictions.rds"))

message("[04b] done.")
