## ---------------------------------------------------------------------------
## 00_setup.R
## Spatio-temporal and machine-learning analysis of under-five mortality
## in Nepal, NDHS 2011 / 2016 / 2022
##
## Global options, paths, package loading and plotting theme.
## Sourced by every downstream script.
## ---------------------------------------------------------------------------

options(stringsAsFactors = FALSE, scipen = 999, dplyr.summarise.inform = FALSE)

suppressPackageStartupMessages({
  library(haven)       
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(scales)
  library(viridis)
  library(patchwork)
  library(survey)      # design-based (weighted) estimation
  library(spdep)       # Moran's I, LISA, Getis-Ord
})

## ---- Reproducibility ------------------------------------------------------
SEED <- 20260101
set.seed(SEED)

## ---- Paths ----------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## Project root: override with U5M_PROJ if the tree is moved.
PROJ <- Sys.getenv("U5M_PROJ",
                   "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal")
if (!dir.exists(PROJ)) PROJ <- getwd()

## Input data (not distributed): override with U5M_DHS_DIR and U5M_SHP_DIR.
DHS_DIR <- Sys.getenv("U5M_DHS_DIR", "D:/khoj_nso/khojibesh")
SHP_DIR <- Sys.getenv("U5M_SHP_DIR",
                      "D:/khoj_nso/Articles/Mortality Estimation/Map for mortality estimation/dist_prov_nep_local")

DIR_DERIVED <- file.path(PROJ, "data_derived")
DIR_RESULTS <- file.path(PROJ, "results")
DIR_FIG     <- file.path(PROJ, "figures")
for (d in c(DIR_DERIVED, DIR_RESULTS, DIR_FIG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## Birth-recode (BR) and GPS (GE) files, by survey round
BR_FILES <- c(
  "2011" = file.path(DHS_DIR, "NP_2011_DHS_05042026_534_208568/NPBR61DT/NPBR61FL.DTA"),
  "2016" = file.path(DHS_DIR, "NP_2016_DHS_05042026_535_208568/NPBR7HDT/NPBR7HFL.DTA"),
  "2022" = file.path(DHS_DIR, "NP_2022_DHS_05042026_536_208568/NPBR82DT/NPBR82FL.DTA")
)
GE_FILES <- c(
  "2011" = file.path(DHS_DIR, "NDHS_2011_GPS/NPGE61FL/NPGE61FL.shp"),
  "2016" = file.path(DHS_DIR, "NDHS_2016_GPS/NPGE7AFL/NPGE7AFL.shp"),
  "2022" = file.path(DHS_DIR, "NDHS_2022_GPS/NPGE82FL/NPGE82FL.shp")
)
SHP_DISTRICT <- file.path(SHP_DIR, "dist_final.shp")
SHP_PROVINCE <- file.path(SHP_DIR, "province.shp")

## ---- Analytic design constants -------------------------------------------
## Retrospective window: live births in the 120 months preceding interview.
WINDOW_MONTHS <- 120

## Standard DHS life-table age segments (months) used for the life table and for
## the person-period (discrete-time survival) expansion.
AGE_SEGMENTS <- data.frame(
  seg    = 1:8,
  start  = c(0, 1,  3,  6, 12, 24, 36, 48),
  end    = c(1, 3,  6, 12, 24, 36, 48, 60),
  label  = c("0", "1-2", "3-5", "6-11", "12-23", "24-35", "36-47", "48-59")
)
AGE_SEGMENTS$width <- AGE_SEGMENTS$end - AGE_SEGMENTS$start

## Calendar birth periods for the temporal dimension
PERIOD_BREAKS <- c(2002, 2007, 2012, 2017, 2023)   # left-closed
PERIOD_LABELS <- c("2002-2006", "2007-2011", "2012-2016", "2017-2022")

## CMC helpers ---------------------------------------------------------------
## IMPORTANT: every century-month code (CMC) in the Nepal DHS recodes -- v008
## (date of interview), v011 (mother's date of birth) and b3 (child's date of
## birth) -- is expressed in the *Bikram Sambat* (BS) calendar, not the
## Gregorian calendar. CMC_BS = (BS year - 1900) * 12 + BS month, with month 1
## = Baisakh, which begins around 13-14 April of the Gregorian year BS - 57.
##
## Differences between two CMCs (ages, exposure, birth intervals) are
## calendar-invariant and need no conversion; only calendar labelling does.
cmc_bs_year  <- function(cmc) 1900 + floor((cmc - 1) / 12)
cmc_bs_month <- function(cmc) ((cmc - 1) %% 12) + 1

## BS CMC -> Gregorian decimal year. The 0.2849 offset is the fraction of the
## Gregorian year elapsed at 14 April (day 104 of 365).
bs_cmc_to_ad <- function(cmc) {
  (cmc_bs_year(cmc) - 57) + 0.2849 + (cmc_bs_month(cmc) - 0.5) / 12
}
bs_cmc_to_ad_year <- function(cmc) floor(bs_cmc_to_ad(cmc))

## ---- Plot theme -----------------------------------------------------------
theme_u5m <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      plot.title       = element_text(face = "bold", size = base_size + 2),
      plot.subtitle    = element_text(colour = "grey30", size = base_size - 1),
      plot.caption     = element_text(colour = "grey40", size = base_size - 2, hjust = 0),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      legend.position  = "bottom",
      legend.key.height = unit(0.4, "cm")
    )
}
theme_map <- function(base_size = 11) {
  theme_void(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(colour = "grey30", size = base_size - 1),
      strip.text    = element_text(face = "bold", size = base_size),
      legend.position = "right"
    )
}
theme_set(theme_u5m())

save_fig <- function(plot, name, width = 9, height = 7, dpi = 400) {
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), plot,
         width = width, height = height, dpi = dpi, bg = "white")
  ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), plot,
         width = width, height = height, bg = "white")
  invisible(plot)
}

message("[00_setup] project = ", PROJ)
