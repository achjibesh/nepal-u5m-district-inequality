## ---------------------------------------------------------------------------
## 03_spatial_analysis.R
##
## Exploratory spatial and spatio-temporal statistics on the 77 x 4
## district-period lattice:
##   * empirical-Bayes smoothing of sparse district rates
##   * global spatial autocorrelation (Moran's I, Geary's C, and the
##     Assuncao-Reis empirical-Bayes index)
##   * local indicators of spatial association (LISA) with FDR control
##   * Getis-Ord Gi* hot- and cold-spots
##   * a Kulldorff space-time scan statistic (discrete Poisson model,
##     cylindrical windows, Monte Carlo inference)
##
## Outputs (results/): moran_global.csv, lisa_clusters.csv,
##                     lisa_persistent_hh.csv, getis_gi.csv,
##                     scan_clusters.csv, scan_membership.csv
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))
suppressPackageStartupMessages(library(Matrix))

district_sf <- readRDS(file.path(DIR_DERIVED, "district_sf.rds")) %>% arrange(dcode)
dp          <- readRDS(file.path(DIR_DERIVED, "district_period.rds")) %>%
  arrange(period, dcode)
stopifnot(all(sort(unique(dp$dcode)) == district_sf$dcode))

## ===========================================================================
## 1. Spatial weights
## ===========================================================================
## Queen contiguity, row-standardised. Nepal's district lattice is fully
## connected, so no island correction is required.

nb <- poly2nb(district_sf, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
message("[03] neighbour links: ", sum(card(nb)),
        " | mean = ", round(mean(card(nb)), 2),
        " | districts with no neighbour: ", sum(card(nb) == 0))
saveRDS(list(nb = nb, listw = lw), file.path(DIR_DERIVED, "spatial_weights.rds"))

## ===========================================================================
## 2. Empirical-Bayes smoothing of district rates
## ===========================================================================
## A district-period cell contains on average only about 100 sampled births
## and a handful of deaths, so the raw life-table rate is dominated by
## sampling noise. Following standard practice for rare-event rates on
## unequal denominators, every autocorrelation statistic is reported both on
## the raw rate and on the Marshall global empirical-Bayes estimate, which
## shrinks each district towards the national rate in inverse proportion to
## the precision of its own denominator.

eb_smooth <- function(deaths, births) {
  ok <- births > 0
  eb <- rep(NA_real_, length(deaths))
  eb[ok] <- spdep::EBest(n = deaths[ok], x = births[ok])$estmm * 1000
  eb[!ok] <- 1000 * sum(deaths) / sum(births)   # unsampled cells: national rate
  eb
}

dp <- dp %>%
  group_by(period) %>%
  mutate(u5mr_eb = eb_smooth(deaths, births)) %>%
  ungroup() %>%
  arrange(period, dcode)
saveRDS(dp, file.path(DIR_DERIVED, "district_period.rds"))

## District rates pooled over 2002-2022 (about 400 sampled births per district)
dp_pool <- dp %>%
  group_by(dcode, district, province) %>%
  summarise(births = sum(births), deaths = sum(deaths),
            expected = sum(expected, na.rm = TRUE), .groups = "drop") %>%
  arrange(dcode) %>%
  mutate(crude   = 1000 * deaths / pmax(births, 1),
         u5mr_eb = eb_smooth(deaths, births),
         smr     = deaths / pmax(expected, 1e-6))
saveRDS(dp_pool, file.path(DIR_DERIVED, "district_pooled.rds"))
write.csv(dp_pool, file.path(DIR_RESULTS, "district_pooled.csv"), row.names = FALSE)

## ===========================================================================
## 3. Global spatial autocorrelation
## ===========================================================================

moran_one <- function(x, label, deaths = NULL, births = NULL) {
  mi <- moran.mc(x, lw, nsim = 9999, zero.policy = TRUE)
  gc <- geary.mc(x, lw, nsim = 9999, zero.policy = TRUE)
  out <- tibble::tibble(stratum = label,
                        moran_I = unname(mi$statistic), moran_p = mi$p.value,
                        geary_C = unname(gc$statistic), geary_p = gc$p.value,
                        eb_index = NA_real_, eb_p = NA_real_)
  if (!is.null(deaths)) {
    ebi <- EBImoran.mc(n = deaths, x = pmax(births, 1), listw = lw,
                       nsim = 9999, zero.policy = TRUE)
    out$eb_index <- unname(ebi$statistic); out$eb_p <- ebi$p.value
  }
  out
}

global_moran <- bind_rows(
  moran_one(dp_pool$u5mr_eb, "2002-2022 (pooled)",
            dp_pool$deaths, dp_pool$births) %>% mutate(rate = "Empirical-Bayes"),
  moran_one(dp_pool$crude, "2002-2022 (pooled)") %>% mutate(rate = "Raw"),
  map_dfr(PERIOD_LABELS, function(p) {
    d <- dp %>% filter(period == p)
    bind_rows(
      moran_one(d$u5mr_eb, p, d$deaths, d$births) %>% mutate(rate = "Empirical-Bayes"),
      moran_one(tidyr::replace_na(d$crude_rate, 0), p) %>% mutate(rate = "Raw")
    )
  })
) %>%
  select(stratum, rate, moran_I, moran_p, geary_C, geary_p, eb_index, eb_p)

print(as.data.frame(global_moran), digits = 3)
write.csv(global_moran, file.path(DIR_RESULTS, "moran_global.csv"), row.names = FALSE)

## ===========================================================================
## 4. LISA and Getis-Ord Gi*
## ===========================================================================
## Local Moran's I with 9,999 conditional permutations on the empirical-Bayes
## rate; p-values are adjusted with the Benjamini-Hochberg false discovery
## rate to control the 77 simultaneous tests.

lisa_one <- function(x, keys, label) {
  lm_res <- localmoran_perm(x, lw, nsim = 9999, zero.policy = TRUE)
  z    <- as.numeric(scale(x))
  lagz <- lag.listw(lw, z, zero.policy = TRUE)
  pval <- lm_res[, ncol(lm_res)]
  padj <- p.adjust(pval, method = "BH")
  quad <- case_when(z > 0 & lagz > 0 ~ "High-High",
                    z < 0 & lagz < 0 ~ "Low-Low",
                    z > 0 & lagz < 0 ~ "High-Low",
                    TRUE             ~ "Low-High")
  ## Cluster membership is reported at the conventional p < 0.05 permutation
  ## threshold, as is standard in the applied LISA literature; the column
  ## `cluster_fdr` additionally shows which districts survive Benjamini-
  ## Hochberg control of the 77 simultaneous tests.
  keys %>% mutate(
    stratum = label, value = x,
    Ii = lm_res[, "Ii"], z_Ii = lm_res[, "Z.Ii"],
    p_value = pval, p_fdr = padj,
    cluster     = if_else(pval > 0.05, "Not significant", quad),
    cluster_fdr = if_else(padj > 0.05, "Not significant", quad))
}

keys_pool <- dp_pool %>% select(dcode, district, province)

lisa_all <- bind_rows(
  lisa_one(dp_pool$u5mr_eb, keys_pool, "2002-2022 (pooled)"),
  map_dfr(PERIOD_LABELS, function(p) {
    d <- dp %>% filter(period == p)
    lisa_one(d$u5mr_eb, d %>% select(dcode, district, province), p)
  })
)
write.csv(lisa_all, file.path(DIR_RESULTS, "lisa_clusters.csv"), row.names = FALSE)
print(lisa_all %>% filter(cluster != "Not significant") %>%
        count(stratum, cluster) %>% as.data.frame())

persistent_hh <- lisa_all %>%
  filter(stratum %in% PERIOD_LABELS, cluster == "High-High") %>%
  count(dcode, district, province, name = "n_periods") %>%
  arrange(desc(n_periods), district)
write.csv(persistent_hh, file.path(DIR_RESULTS, "lisa_persistent_hh.csv"), row.names = FALSE)
print(as.data.frame(persistent_hh))

nb_star <- include.self(nb)
lw_star <- nb2listw(nb_star, style = "B", zero.policy = TRUE)

getis_one <- function(x, keys, label) {
  gz <- as.numeric(localG(x, lw_star, zero.policy = TRUE))
  keys %>%
    mutate(stratum = label, value = x, gi_star = gz,
           p_two = 2 * pnorm(abs(gz), lower.tail = FALSE)) %>%
    mutate(p_fdr = p.adjust(p_two, "BH"),
           hotspot = case_when(p_two > 0.05 ~ "Not significant",
                               gi_star > 0  ~ "Hot spot",
                               TRUE         ~ "Cold spot"),
           hotspot_fdr = case_when(p_fdr > 0.05 ~ "Not significant",
                                   gi_star > 0  ~ "Hot spot",
                                   TRUE         ~ "Cold spot"))
}

getis <- bind_rows(
  getis_one(dp_pool$u5mr_eb, keys_pool, "2002-2022 (pooled)"),
  map_dfr(PERIOD_LABELS, function(p) {
    d <- dp %>% filter(period == p)
    getis_one(d$u5mr_eb, d %>% select(dcode, district, province), p)
  })
)
write.csv(getis, file.path(DIR_RESULTS, "getis_gi.csv"), row.names = FALSE)
print(getis %>% filter(hotspot != "Not significant") %>%
        count(stratum, hotspot) %>% as.data.frame())

## ===========================================================================
## 5. Kulldorff space-time scan statistic (discrete Poisson)
## ===========================================================================
## Cylindrical scanning windows are the Cartesian product of
##   * circular spatial bases: for each district centroid, the k nearest
##     districts for k = 1 .. K, K capped so that a window never covers more
##     than 50% of the nationally expected deaths, and
##   * contiguous temporal intervals covering at most half the study period.
## For each window the discrete Poisson log-likelihood ratio is
##   LLR = c log(c/n) + (C - c) log((C - c)/(C - n))   when c > n, else 0,
## with c, n the observed and expected deaths inside the window and C, N the
## national totals. Expected counts come from indirect standardisation of the
## person-period file by the national age-segment hazard schedule. Inference
## uses 999 Monte Carlo replicates in which the C observed deaths are
## re-allocated multinomially with probabilities proportional to expectation.

centroids <- suppressWarnings(st_coordinates(st_centroid(st_geometry(district_sf))))
D   <- as.matrix(dist(centroids))
ord <- t(apply(D, 1, order))                       # nearest-first ordering

cell <- dp %>%
  mutate(expected = tidyr::replace_na(expected, 0)) %>%
  arrange(period_num, dcode)

n_dist <- nrow(district_sf)
n_per  <- length(PERIOD_LABELS)
cell_index <- matrix(seq_len(n_dist * n_per), nrow = n_dist, ncol = n_per)

obs_vec <- cell$deaths
## Kulldorff's discrete Poisson likelihood assumes that the expected counts
## sum to the observed total, so the indirectly standardised expectations are
## rescaled by C / sum(E) before scanning.
exp_vec <- cell$expected * sum(cell$deaths) / sum(cell$expected)
C_tot   <- sum(obs_vec)
N_tot   <- sum(exp_vec)
MAXFRAC <- 0.5
message("[03] total observed deaths = ", C_tot,
        " | total expected (rescaled) = ", round(N_tot, 1))

time_intervals <- do.call(rbind, lapply(1:n_per, function(a)
  do.call(rbind, lapply(a:n_per, function(b) c(a, b)))))
time_intervals <- time_intervals[
  (time_intervals[, 2] - time_intervals[, 1] + 1) <= ceiling(n_per * MAXFRAC), ,
  drop = FALSE]

rows <- integer(0); cols <- integer(0)
meta <- list(); cyl_cells <- list(); cyl <- 0L

for (i in seq_len(n_dist)) {
  members <- integer(0)
  for (k in seq_len(n_dist)) {
    members <- c(members, ord[i, k])
    if (sum(exp_vec[as.vector(cell_index[members, ])]) > MAXFRAC * N_tot) break
    for (t in seq_len(nrow(time_intervals))) {
      t1 <- time_intervals[t, 1]; t2 <- time_intervals[t, 2]
      idx <- as.vector(cell_index[members, t1:t2])
      if (sum(exp_vec[idx]) > MAXFRAC * N_tot) next
      cyl <- cyl + 1L
      rows <- c(rows, rep.int(cyl, length(idx)))
      cols <- c(cols, idx)
      cyl_cells[[cyl]] <- idx
      meta[[cyl]] <- c(centre = i, k = k, t1 = t1, t2 = t2)
    }
  }
}
W    <- sparseMatrix(i = rows, j = cols, x = 1, dims = c(cyl, n_dist * n_per))
meta <- as.data.frame(do.call(rbind, meta))
message("[03] scan cylinders evaluated: ", nrow(meta))

llr_fun <- function(cvec, nvec, C = C_tot, N = N_tot) {
  out <- numeric(length(cvec))
  ok  <- is.finite(cvec) & is.finite(nvec) &
    cvec > nvec & cvec > 0 & nvec > 0 & (C - cvec) > 0 & (N - nvec) > 0
  cc  <- cvec[ok]; nn <- nvec[ok]
  val <- cc * log(cc / nn) + (C - cc) * log((C - cc) / (N - nn))
  val[!is.finite(val)] <- 0
  out[ok] <- val
  out
}

n_win <- as.vector(W %*% exp_vec)
c_win <- as.vector(W %*% obs_vec)
llr   <- llr_fun(c_win, n_win)

set.seed(SEED)
NSIM   <- 999
p_cell <- exp_vec / N_tot
sim    <- matrix(rmultinom(NSIM, size = C_tot, prob = p_cell), ncol = NSIM)
sim_c  <- as.matrix(W %*% sim)
max_llr_sim <- apply(sim_c, 2, function(cv) {
  v <- llr_fun(cv, n_win)
  if (!length(v)) 0 else max(v, 0, na.rm = TRUE)
})
message("[03] observed max LLR = ", round(max(llr), 2),
        " | Monte Carlo max LLR: median ", round(median(max_llr_sim), 2),
        ", 95th pct ", round(quantile(max_llr_sim, 0.95), 2),
        " | non-finite replicates: ", sum(!is.finite(max_llr_sim)))

## Most likely and secondary clusters, constrained not to overlap in space.
ordl <- order(llr, decreasing = TRUE)
selected <- integer(0); used <- integer(0)
for (j in ordl) {
  if (llr[j] <= 0) break
  dists_j <- unique(((cyl_cells[[j]] - 1L) %% n_dist) + 1L)
  if (length(intersect(dists_j, used)) > 0) next
  selected <- c(selected, j)
  used <- union(used, dists_j)
  if (length(selected) >= 8) break
}

cluster_row <- function(r) {
  j <- selected[r]
  dists_j <- sort(unique(((cyl_cells[[j]] - 1L) %% n_dist) + 1L))
  ## values are pulled out before the tibble is built: a column named `llr`
  ## would otherwise mask the global vector for the columns after it
  llr_j <- llr[j]
  p_j   <- (1 + sum(max_llr_sim >= llr_j)) / (NSIM + 1)
  rr_j  <- (c_win[j] / n_win[j]) / ((C_tot - c_win[j]) / (N_tot - n_win[j]))
  tibble::tibble(
    rank        = r,
    type        = if (r == 1) "Most likely" else "Secondary",
    period_from = PERIOD_LABELS[meta$t1[j]],
    period_to   = PERIOD_LABELS[meta$t2[j]],
    n_districts = length(dists_j),
    districts   = paste(sort(district_sf$district[dists_j]), collapse = ", "),
    provinces   = paste(sort(unique(as.character(district_sf$province[dists_j]))),
                        collapse = ", "),
    observed = c_win[j],
    expected = round(n_win[j], 1),
    rr       = round(rr_j, 2),
    llr      = round(llr_j, 2),
    p_value  = p_j)
}

scan_clusters <- map_dfr(seq_along(selected), cluster_row)
print(scan_clusters %>%
        select(rank, type, period_from, period_to, n_districts, observed,
               expected, rr, llr, p_value) %>% as.data.frame())
write.csv(scan_clusters, file.path(DIR_RESULTS, "scan_clusters.csv"), row.names = FALSE)

scan_membership <- map_dfr(seq_along(selected), function(r) {
  j <- selected[r]
  dists_j <- sort(unique(((cyl_cells[[j]] - 1L) %% n_dist) + 1L))
  rr_j <- (c_win[j] / n_win[j]) / ((C_tot - c_win[j]) / (N_tot - n_win[j]))
  p_j  <- (1 + sum(max_llr_sim >= llr[j])) / (NSIM + 1)
  tibble::tibble(cluster_rank = r,
                 dcode = district_sf$dcode[dists_j],
                 district = district_sf$district[dists_j],
                 period_from = PERIOD_LABELS[meta$t1[j]],
                 period_to   = PERIOD_LABELS[meta$t2[j]],
                 rr = round(rr_j, 2),
                 p_value = p_j)
})
saveRDS(scan_membership, file.path(DIR_DERIVED, "scan_membership.rds"))
write.csv(scan_membership, file.path(DIR_RESULTS, "scan_membership.csv"), row.names = FALSE)

message("[03] done.")
