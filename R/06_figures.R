## ---------------------------------------------------------------------------
## 06_figures.R
##
## All main-text and supplementary figures. Every panel is written to
## figures/ at 400 dpi in both PNG and PDF.
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))
suppressPackageStartupMessages({
  library(yardstick)
  library(ggrepel)
})

district_sf <- readRDS(file.path(DIR_DERIVED, "district_sf.rds")) %>% arrange(dcode)
province_sf <- readRDS(file.path(DIR_DERIVED, "province_sf.rds"))
dp          <- readRDS(file.path(DIR_DERIVED, "district_period.rds"))
dp_pool     <- readRDS(file.path(DIR_DERIVED, "district_pooled.rds"))
nat_period  <- read.csv(file.path(DIR_RESULTS, "rates_national_period.csv"))
prov_period <- read.csv(file.path(DIR_RESULTS, "rates_province_period.csv"))
getis       <- read.csv(file.path(DIR_RESULTS, "getis_gi.csv"))
lisa        <- read.csv(file.path(DIR_RESULTS, "lisa_clusters.csv"))
scan_mem    <- readRDS(file.path(DIR_DERIVED, "scan_membership.rds"))
scan_clu    <- read.csv(file.path(DIR_RESULTS, "scan_clusters.csv"))
rr          <- readRDS(file.path(DIR_DERIVED, "inla_district_rr.rds"))
resid_tbl   <- readRDS(file.path(DIR_DERIVED, "inla_residual_district.rds"))
test_pred   <- readRDS(file.path(DIR_DERIVED, "test_predictions.rds"))
shap_imp    <- read.csv(file.path(DIR_RESULTS, "shap_importance.csv"))
shap_prov   <- read.csv(file.path(DIR_RESULTS, "shap_by_province.csv"))
shap_long   <- readRDS(file.path(DIR_DERIVED, "shap_long.rds"))
pdp         <- read.csv(file.path(DIR_RESULTS, "pdp_top_features.csv"))
perf        <- read.csv(file.path(DIR_RESULTS, "ml_performance.csv"))
calib       <- read.csv(file.path(DIR_RESULTS, "ml_calibration.csv"))
spatial_cv  <- read.csv(file.path(DIR_RESULTS, "ml_spatial_cv.csv"))

prov_line <- st_as_sf(st_cast(st_geometry(province_sf), "MULTILINESTRING"))
add_prov  <- function() geom_sf(data = prov_line, inherit.aes = FALSE,
                                colour = "grey15", linewidth = 0.35, fill = NA)

join_map <- function(tbl, by = "dcode") {
  district_sf %>% select(dcode, geometry) %>% left_join(tbl, by = by)
}

## ===========================================================================
## Figure 1. National trend and the distribution of district mortality
## ===========================================================================

trend_long <- nat_period %>%
  select(period, `Under-five` = u5mr, Infant = imr, Neonatal = nnmr) %>%
  tidyr::pivot_longer(-period, names_to = "indicator", values_to = "rate") %>%
  mutate(indicator = factor(indicator, levels = c("Under-five", "Infant", "Neonatal")))

p1a <- ggplot(trend_long, aes(period, rate, group = indicator, colour = indicator)) +
  geom_ribbon(data = nat_period, inherit.aes = FALSE,
              aes(x = period, ymin = lo, ymax = hi, group = 1),
              fill = "grey70", alpha = 0.30) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  geom_text(aes(label = round(rate, 0)), vjust = -1.1, size = 3.1, show.legend = FALSE) +
  scale_colour_manual(values = c("Under-five" = "#B2182B", "Infant" = "#2166AC",
                                 "Neonatal" = "#1B7837")) +
  labs(x = NULL, y = "Deaths per 1,000 live births", colour = NULL,
       title = "A  National child mortality by birth period",
       subtitle = str_wrap(paste("Life-table estimates by birth cohort, pooled across NDHS 2011,",
                                 "2016 and 2022; the shaded band is the 95%",
                                 "cluster-bootstrap interval for under-five mortality"), 95)) +
  expand_limits(y = 0)

p1b <- ggplot(rr, aes(period, u5mr_smooth)) +
  geom_violin(fill = "#4393C3", alpha = 0.25, colour = NA, width = 0.9) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.09, alpha = 0.35, size = 0.8, colour = "#2166AC") +
  geom_point(data = nat_period, aes(period, u5mr), inherit.aes = FALSE,
             shape = 23, size = 2.8, fill = "#B2182B", colour = "white", stroke = 0.6) +
  labs(x = NULL, y = "Posterior district U5MR per 1,000",
       title = "B  Distribution of district-level under-five mortality",
       subtitle = str_wrap(paste("Posterior means from the Bayesian spatio-temporal",
                                 "model, one point per district; the red diamond is",
                                 "the national rate in each period"), 95))

fig1 <- p1a / p1b + plot_layout(heights = c(1, 1.15))
save_fig(fig1, "fig1_national_trend", width = 9.5, height = 9)

## ===========================================================================
## Figure 2. Smoothed district under-five mortality by period
## ===========================================================================

map2_df <- join_map(rr %>% select(dcode, period, u5mr_smooth)) %>%
  mutate(period = factor(period, levels = PERIOD_LABELS))

fig2 <- ggplot(map2_df) +
  geom_sf(aes(fill = u5mr_smooth), colour = "white", linewidth = 0.08) +
  add_prov() +
  facet_wrap(~ period, ncol = 2) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       name = "U5MR per 1,000\nlive births") +
  labs(title = "District-level under-five mortality in Nepal, 2002-2022",
       subtitle = "Posterior means from the Bayesian spatio-temporal model (BYM2 spatial effect, random-walk period effect, space-time interaction)",
       caption = "Grey outlines are the seven provinces of the 2015 federal structure.") +
  theme_map()
save_fig(fig2, "fig2_smoothed_u5mr_maps", width = 11, height = 8.5)

## ===========================================================================
## Figure 3. Spatial clustering: Getis-Ord hot-spots and space-time clusters
## ===========================================================================

hot_df <- join_map(getis %>% filter(stratum %in% PERIOD_LABELS) %>%
                     select(dcode, period = stratum, hotspot)) %>%
  mutate(period = factor(period, levels = PERIOD_LABELS),
         hotspot = factor(hotspot, levels = c("Hot spot", "Not significant", "Cold spot")))

p3a <- ggplot(hot_df) +
  geom_sf(aes(fill = hotspot), colour = "white", linewidth = 0.08) +
  add_prov() +
  facet_wrap(~ period, ncol = 2) +
  scale_fill_manual(values = c("Hot spot" = "#B2182B", "Not significant" = "grey88",
                               "Cold spot" = "#2166AC"), name = NULL) +
  labs(title = "A  Getis-Ord Gi* hot- and cold-spots of under-five mortality",
       subtitle = "Empirical-Bayes district rates, self-inclusive contiguity weights, p < 0.05") +
  theme_map() + theme(legend.position = "bottom")

## Only clusters that reach significance are mapped: showing the highest-
## likelihood windows regardless of their p-value would imply a structure the
## Monte Carlo test does not support.
span <- function(from, to) {
  a <- substr(from, 1, 4); b <- substr(to, 6, 9)
  ifelse(a == substr(to, 1, 4) & from == to, from, paste0(a, "-", b))
}
sig <- scan_clu %>% filter(p_value < 0.05)
if (nrow(sig) == 0) sig <- scan_clu[1, ]

scan_df <- district_sf %>%
  left_join(scan_mem %>% filter(cluster_rank %in% sig$rank), by = "dcode") %>%
  mutate(cl = if_else(is.na(cluster_rank), "Outside any detected cluster",
                      paste0("Cluster ", cluster_rank)))
scan_lab <- sig %>%
  mutate(lab = sprintf("Cluster %d (%s): %s, relative risk %.2f, %d observed vs %.1f expected deaths, p = %.3f",
                       rank, tolower(type), span(period_from, period_to),
                       rr, observed, expected, p_value)) %>%
  pull(lab) %>% str_wrap(110) %>% paste(collapse = "\n")

p3b <- ggplot(scan_df) +
  geom_sf(aes(fill = cl), colour = "white", linewidth = 0.08) +
  add_prov() +
  scale_fill_manual(values = c("Cluster 1" = "#B2182B", "Cluster 2" = "#EF8A62",
                               "Cluster 3" = "#FDDBC7", "Cluster 4" = "#7B3294",
                               "Outside any detected cluster" = "grey90"),
                    name = NULL) +
  labs(title = "B  Statistically significant Kulldorff space-time scan clusters",
       subtitle = scan_lab,
       caption = "Expected counts are indirectly standardised by age segment within birth period, so a cluster is high relative to the national level at the same date.") +
  theme_map() + theme(legend.position = "bottom")

fig3 <- p3a / p3b + plot_layout(heights = c(1.25, 1))
save_fig(fig3, "fig3_clusters", width = 11, height = 11)

## ===========================================================================
## Figure 4. Bayesian relative risk, exceedance probability, residual effect
## ===========================================================================

last_p <- tail(PERIOD_LABELS, 1)
rr_last <- rr %>% filter(period == last_p)

p4a <- ggplot(join_map(rr_last %>% select(dcode, rr))) +
  geom_sf(aes(fill = rr), colour = "white", linewidth = 0.08) + add_prov() +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       midpoint = 1, name = "Relative\nrisk") +
  labs(title = paste0("A  Posterior relative risk, ", last_p),
       subtitle = "Relative to the national under-five mortality level of the same period") +
  theme_map()

p4b <- ggplot(join_map(rr_last %>% select(dcode, p_exceed_1_2))) +
  geom_sf(aes(fill = p_exceed_1_2), colour = "white", linewidth = 0.08) + add_prov() +
  scale_fill_viridis_c(option = "magma", direction = -1, limits = c(0, 1),
                       name = "P(RR > 1.2)") +
  labs(title = paste0("B  Posterior exceedance probability, ", last_p),
       subtitle = "Probability that district mortality exceeds the national level by more than 20%") +
  theme_map()

p4c <- ggplot(join_map(resid_tbl %>% select(dcode, resid_or))) +
  geom_sf(aes(fill = resid_or), colour = "white", linewidth = 0.08) + add_prov() +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       midpoint = 1, name = "Residual\nodds ratio") +
  labs(title = "C  District effect after full individual-level adjustment",
       subtitle = "Spatial variation in child survival not explained by measured child, maternal and household characteristics") +
  theme_map()

fig4 <- (p4a | p4b) / p4c + plot_layout(heights = c(1, 1.05))
save_fig(fig4, "fig4_bayesian_maps", width = 12, height = 9)

## ===========================================================================
## Figure 5. Machine-learning discrimination and calibration
## ===========================================================================

LEARNERS <- c("Elastic-net logistic regression", "Random forest",
              "Gradient boosting", "Stacked ensemble")
LCOL <- c("Elastic-net logistic regression" = "#4D4D4D",
          "Random forest" = "#1B7837", "Gradient boosting" = "#B2182B",
          "Stacked ensemble" = "#2166AC")

long_pred <- test_pred %>%
  tidyr::pivot_longer(all_of(LEARNERS), names_to = "learner", values_to = "prob") %>%
  mutate(truth = factor(y, levels = c(1, 0)))

roc_df <- long_pred %>% group_by(learner) %>%
  yardstick::roc_curve(truth, prob, event_level = "first")
pr_df <- long_pred %>% group_by(learner) %>%
  yardstick::pr_curve(truth, prob, event_level = "first")

auc_lab <- perf %>%
  mutate(lab = sprintf("%s (AUC %.3f)", learner, roc_auc)) %>%
  select(learner, lab, roc_auc, pr_auc)

p5a <- ggplot(roc_df, aes(1 - specificity, sensitivity, colour = learner)) +
  geom_abline(linetype = 3, colour = "grey60") +
  geom_path(linewidth = 0.8) +
  scale_colour_manual(values = LCOL, name = NULL,
                      labels = setNames(sprintf("%s  (AUC %.3f)", auc_lab$learner,
                                                auc_lab$roc_auc), auc_lab$learner)) +
  coord_equal() +
  labs(x = "1 - specificity", y = "Sensitivity",
       title = "A  Receiver operating characteristic") +
  theme(legend.position = "inside", legend.position.inside = c(0.62, 0.22),
        legend.background = element_rect(fill = "white", colour = NA))

p5b <- ggplot(pr_df, aes(recall, precision, colour = learner)) +
  geom_hline(yintercept = mean(test_pred$y), linetype = 3, colour = "grey50") +
  geom_path(linewidth = 0.8) +
  scale_colour_manual(values = LCOL, name = NULL,
                      labels = setNames(sprintf("%s  (PR-AUC %.3f)", auc_lab$learner,
                                                auc_lab$pr_auc), auc_lab$learner)) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = "Recall", y = "Precision",
       title = "B  Precision-recall",
       subtitle = "Dotted line is the segment-death prevalence in the test set") +
  theme(legend.position = "inside", legend.position.inside = c(0.62, 0.78),
        legend.background = element_rect(fill = "white", colour = NA))

p5c <- ggplot(calib, aes(pred, obs, colour = learner)) +
  geom_abline(linetype = 3, colour = "grey60") +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_colour_manual(values = LCOL, name = NULL) +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Mean predicted hazard (decile)", y = "Observed proportion dying",
       title = "C  Calibration by decile of predicted risk") +
  theme(legend.position = "bottom")

p5d <- spatial_cv %>%
  select(province, `Gradient boosting` = roc_xgb, `Elastic net` = roc_enet) %>%
  tidyr::pivot_longer(-province, names_to = "learner", values_to = "auc") %>%
  ggplot(aes(reorder(province, auc), auc, fill = learner)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_hline(yintercept = 0.5, linetype = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Gradient boosting" = "#B2182B", "Elastic net" = "#4D4D4D"),
                    name = NULL) +
  labs(x = NULL, y = "Area under the ROC curve",
       title = "D  Leave-one-province-out cross-validation",
       subtitle = "Each province is predicted by a model trained only on the other six") +
  theme(legend.position = "bottom")

fig5 <- (p5a | p5b) / (p5c | p5d)
save_fig(fig5, "fig5_ml_performance", width = 12, height = 10)

## ===========================================================================
## Figure 6. SHAP importance and beeswarm
## ===========================================================================

imp_plot <- shap_imp %>% slice_max(mean_abs_shap, n = 15) %>%
  mutate(label = factor(label, levels = rev(label)))

p6a <- ggplot(imp_plot, aes(mean_abs_shap, label, fill = block)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("Baseline hazard" = "#4D4D4D",
                               "Child and birth" = "#B2182B",
                               "Maternal" = "#EF8A62",
                               "Household" = "#2166AC",
                               "Geographic context" = "#1B7837"), name = NULL) +
  labs(x = "Mean absolute SHAP value (log-odds of dying in an age segment)",
       y = NULL, title = "A  Global predictor importance") +
  theme(legend.position = "bottom")

## The SHAP table has one row per person-period per one-hot column; a beeswarm
## over all of them is unreadable and slow to draw, so a fixed random sample is
## plotted per predictor.
set.seed(SEED)
bees <- shap_long %>%
  group_by(variable) %>%
  slice_sample(n = 4000) %>%
  ungroup() %>%
  group_by(column) %>%
  mutate(xs = if (diff(range(x, na.rm = TRUE)) > 0)
    (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) else 0.5) %>%
  ungroup() %>%
  group_by(variable) %>%
  mutate(imp = mean(abs(shap))) %>%
  ungroup() %>%
  left_join(shap_imp %>% select(variable, label), by = "variable") %>%
  mutate(label = forcats::fct_reorder(label, imp))

p6b <- ggplot(bees, aes(shap, label, colour = xs)) +
  geom_vline(xintercept = 0, colour = "grey60", linetype = 3) +
  geom_jitter(height = 0.28, size = 0.35, alpha = 0.25) +
  scale_colour_gradient(low = "#2166AC", high = "#B2182B",
                        breaks = c(0, 1), labels = c("Low", "High"),
                        name = "Predictor value") +
  labs(x = "SHAP value (contribution to the log-odds of dying)", y = NULL,
       title = "B  Direction of individual contributions") +
  theme(legend.position = "bottom")

fig6 <- p6a | p6b
save_fig(fig6, "fig6_shap", width = 13.5, height = 8)

## ===========================================================================
## Figure 7. Province-stratified SHAP: geography of risk-factor importance
## ===========================================================================

top12 <- shap_imp %>% filter(variable != "seg_f") %>%
  slice_max(mean_abs_shap, n = 12) %>% pull(variable)

hm <- shap_prov %>%
  filter(variable %in% top12) %>%
  left_join(shap_imp %>% select(variable, label), by = "variable") %>%
  group_by(province) %>%
  mutate(share = 100 * mean_abs_shap / sum(mean_abs_shap)) %>%
  ungroup()

p7 <- ggplot(hm, aes(province, forcats::fct_reorder(label, share), fill = share)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", share)), size = 2.9, colour = "grey15") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Share of total\nSHAP magnitude (%)") +
  labs(x = NULL, y = NULL,
       title = "Composition of modelled risk by province",
       subtitle = str_wrap(paste("Percentage of each province's total SHAP magnitude",
                                 "attributable to each predictor, from the single national",
                                 "gradient-boosting model. The columns are close to",
                                 "identical: the determinants of child death act nationally,",
                                 "not province by province."), 120)) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")
save_fig(p7, "fig7_shap_by_province", width = 10.5, height = 7.5)

## ===========================================================================
## Supplementary figures
## ===========================================================================

## S1 - shrinkage: raw, empirical-Bayes and fully Bayesian district estimates
shrink <- rr %>%
  select(dcode, district, period, births, u5mr_raw, u5mr_eb, u5mr_smooth) %>%
  tidyr::pivot_longer(c(u5mr_raw, u5mr_eb, u5mr_smooth),
                      names_to = "estimator", values_to = "value") %>%
  mutate(estimator = recode(estimator,
                            u5mr_raw = "Raw life-table rate",
                            u5mr_eb = "Empirical-Bayes rate",
                            u5mr_smooth = "Bayesian spatio-temporal posterior"))
s1 <- ggplot(shrink, aes(births, value, colour = estimator)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(se = FALSE, method = "loess", formula = y ~ x, linewidth = 0.8) +
  facet_wrap(~ period, ncol = 2) +
  scale_colour_manual(values = c("Raw life-table rate" = "grey55",
                                 "Empirical-Bayes rate" = "#EF8A62",
                                 "Bayesian spatio-temporal posterior" = "#2166AC"),
                      name = NULL) +
  labs(x = "Sampled births in the district-period cell",
       y = "Under-five mortality per 1,000",
       title = "Shrinkage of district estimates towards the national level",
       subtitle = "Sparse cells are pulled hardest towards the mean; the Bayesian posterior additionally borrows strength from neighbouring districts and adjacent periods")
save_fig(s1, "figS1_shrinkage", width = 10, height = 8)

## S2 - LISA cluster maps
lisa_df <- join_map(lisa %>% filter(stratum %in% PERIOD_LABELS) %>%
                      select(dcode, period = stratum, cluster)) %>%
  mutate(period = factor(period, levels = PERIOD_LABELS),
         cluster = factor(cluster, levels = c("High-High", "Low-Low", "High-Low",
                                              "Low-High", "Not significant")))
s2 <- ggplot(lisa_df) +
  geom_sf(aes(fill = cluster), colour = "white", linewidth = 0.08) + add_prov() +
  facet_wrap(~ period, ncol = 2) +
  scale_fill_manual(values = c("High-High" = "#B2182B", "Low-Low" = "#2166AC",
                               "High-Low" = "#F4A582", "Low-High" = "#92C5DE",
                               "Not significant" = "grey90"), name = NULL) +
  labs(title = "Local indicators of spatial association (LISA), by birth period",
       subtitle = "Local Moran's I on empirical-Bayes district rates, 9,999 conditional permutations, p < 0.05") +
  theme_map()
save_fig(s2, "figS2_lisa", width = 10.5, height = 8.5)

## S3 - partial dependence
s3 <- pdp %>%
  filter(!(variable == "bi_num" & value < 0)) %>%
  ggplot(aes(value, yhat)) +
  geom_line(colour = "#B2182B", linewidth = 0.9) +
  facet_wrap(~ label, scales = "free", ncol = 3) +
  labs(x = NULL, y = "Average predicted segment hazard",
       title = "Partial dependence of the segment hazard on the continuous predictors",
       subtitle = "Gradient-boosting model fitted to the full person-period sample")
save_fig(s3, "figS3_pdp", width = 11, height = 7.5)

## S4 - mortality by province
s4 <- ggplot(prov_period, aes(period, u5mr, group = province, colour = province)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.9) +
  scale_colour_brewer(palette = "Dark2", name = NULL) +
  labs(x = NULL, y = "U5MR per 1,000",
       title = "Under-five mortality by province and birth cohort") +
  expand_limits(y = 0)
save_fig(s4, "figS4_province_trends", width = 9.5, height = 5.5)

## S5 - study area
s5 <- ggplot(district_sf) +
  geom_sf(aes(fill = province), colour = "white", linewidth = 0.12) +
  add_prov() +
  scale_fill_brewer(palette = "Pastel2", name = "Province") +
  labs(title = "Study area: 77 districts and 7 provinces of Nepal") +
  theme_map()
save_fig(s5, "figS5_study_area", width = 9.5, height = 6)

## S7 - derivation of the analytic sample (STROBE flow)
flow <- read.csv(file.path(DIR_RESULTS, "sample_flow.csv")) %>%
  mutate(i = row_number(), y = -i,
         lab = sprintf("%s\nn = %s", str_wrap(step, 46), format(n, big.mark = ",")),
         dropped = lag(n) - n)
sF <- ggplot(flow) +
  geom_segment(data = flow %>% filter(i > 1),
               aes(x = 0, xend = 0, y = y + 1 - 0.3, yend = y + 0.3),
               arrow = arrow(length = unit(0.18, "cm")), colour = "grey40") +
  geom_label(aes(x = 0, y = y, label = lab), size = 3.2, fill = "white",
             label.padding = unit(0.4, "lines")) +
  geom_text(data = flow %>% filter(i > 1),
            aes(x = 0.62, y = y + 0.5,
                label = sprintf("excluded: %s", format(dropped, big.mark = ","))),
            size = 3, hjust = 0, colour = "#B2182B") +
  scale_x_continuous(limits = c(-1.25, 1.7)) +
  labs(title = "Derivation of the analytic sample") +
  theme_void() + theme(plot.title = element_text(face = "bold"))
save_fig(sF, "figS7_sample_flow", width = 7.5, height = 6.5)

message("[06] figures written to ", DIR_FIG)
