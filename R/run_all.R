## ---------------------------------------------------------------------------
## run_all.R
##
## Runs the whole analysis from the raw NDHS files to the finished manuscript.
##
## Usage, from the project root:
##   Rscript R/run_all.R
##
## Set the environment variable U5M_PROJ if the project tree has been moved.
## Expected wall-clock time on four threads: roughly two to four hours, almost
## all of it in 04_machine_learning.R (hyperparameter search and out-of-fold
## prediction over 209,000 person-periods).
## ---------------------------------------------------------------------------

PROJ <- Sys.getenv("U5M_PROJ",
                   "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal")
Sys.setenv(U5M_PROJ = PROJ)
setwd(PROJ)

STEPS <- c(
  "01_prepare_data.R",          # birth histories -> person-period file
  "02_mortality_rates.R",       # synthetic-cohort life tables
  "03_spatial_analysis.R",      # Moran, LISA, Gi*, space-time scan
  "04_machine_learning.R",      # four learners, validation, SHAP
  "04b_recalibration.R",        # Platt scaling of the class-weighted learners
  "05_bayesian_spatiotemporal.R", # INLA area and individual models
  "06_figures.R",               # all figures
  "06b_figure_exposure.R",      # methods figure: exposure and events
  "07_tables.R"                 # all display tables
)

dir.create(file.path(PROJ, "results"), showWarnings = FALSE, recursive = TRUE)

for (s in STEPS) {
  message("\n=================================================================")
  message("  ", s, "   (", format(Sys.time(), "%H:%M:%S"), ")")
  message("=================================================================")
  t0 <- Sys.time()
  log_file <- file.path(PROJ, "results",
                        paste0("log_", sub("_.*", "", s), ".txt"))
  status <- system2("Rscript", shQuote(file.path(PROJ, "R", s)),
                    stdout = log_file, stderr = log_file)
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
  if (status != 0) {
    stop("Step failed: ", s, " (exit ", status, "). See ", log_file)
  }
  message("  done in ", el, " min")
}

## --- render the manuscript and the supplementary file ----------------------
quarto <- Sys.which("quarto")
if (quarto == "") {
  cand <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
  if (file.exists(cand)) quarto <- cand
}
if (quarto != "") {
  ## renders both documents and assembles manuscript/submission/ with checks
  message("building the submission package")
  system2(file.path(R.home("bin"), "Rscript"),
          shQuote(file.path(PROJ, "manuscript", "build_submission.R")))
} else {
  message("quarto not found on the path; render manuscript/manuscript.qmd manually")
}

message("\nAll steps complete. Outputs in results/, figures/ and manuscript/.")
