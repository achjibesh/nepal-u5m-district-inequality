# District-level inequality in under-five mortality in Nepal, 2002–2022

Analysis code for the manuscript *"District-level inequality in under-five
mortality in Nepal, 2002–2022: a Bayesian spatio-temporal analysis"*
(Jibesh Acharya; submitted to *PLOS ONE*).

**Question.** Where in Nepal are children still dying before their fifth
birthday, has the gap between districts narrowed between 2002 and 2022, and
what determines the risk?

**Data.** Nepal Demographic and Health Surveys 2011, 2016 and 2022 (birth
recode and GPS files), plus a 77-district administrative boundary layer.
Neither is included in this repository (see *Data* below).

---

## Layout

```
├── R/
│   ├── 00_setup.R                  paths, constants, calendar helpers, theme
│   ├── 01_prepare_data.R           birth histories -> person-period file
│   ├── 02_mortality_rates.R        life tables by birth cohort
│   ├── 03_spatial_analysis.R       Moran / LISA / Gi* / space-time scan
│   ├── 04_machine_learning.R       four learners, validation, SHAP
│   ├── 04b_recalibration.R         Platt scaling of the class-weighted models
│   ├── 05_bayesian_spatiotemporal.R  INLA area and individual models
│   ├── 06_figures.R                all figures
│   ├── 06b_figure_exposure.R       exposure figure
│   ├── 07_tables.R                 all display tables
│   └── run_all.R                   driver
├── data_reference/                 district ecological-belt lookup
├── results/                        every statistic as CSV, plus run logs
│   └── tables/                     publication-ready display tables
├── figures/                        PNG and PDF
├── manuscript/
│   ├── manuscript.qmd              main text; pulls numbers from results/
│   ├── additional_file_1_supplementary.qmd   S1 File
│   ├── additional_file_2_STROBE.md           S1 Checklist
│   ├── additional_file_3_TRIPOD_AI.md        S2 Checklist
│   ├── references.bib, plos.csl, reference.docx
│   ├── make_plos_figures.R         PLOS-format TIFFs
│   └── build_plos_submission.R     renders the submission package
└── review/repair_docx.R            Word schema fixes used by the build
```

`data_derived/` (intermediate person-level files) is created by the pipeline
and is deliberately not tracked.

## Data

The NDHS microdata and GPS files cannot be redistributed under the DHS
Program's data-use terms. Register at <https://dhsprogram.com/data/> and
request the Nepal 2011, 2016 and 2022 surveys. The files used are:

| Survey | Birth recode | GPS |
|---|---|---|
| 2011 | `NPBR61FL.DTA` | `NPGE61FL.shp` |
| 2016 | `NPBR7HFL.DTA` | `NPGE7AFL.shp` |
| 2022 | `NPBR82FL.DTA` | `NPGE82FL.shp` |

District and province boundaries (77 districts, 7 provinces) were obtained
from the Survey Department, Government of Nepal (<https://dos.gov.np>).

The expected sub-folder layout under the DHS directory is given in
`BR_FILES` and `GE_FILES` in `R/00_setup.R`; adjust those two vectors if your
download is arranged differently.

## Running it

```bash
export U5M_PROJ=/path/to/this/repo
export U5M_DHS_DIR=/path/to/dhs/downloads
export U5M_SHP_DIR=/path/to/boundaries      # dist_final.shp, province.shp
Rscript R/run_all.R
```

Then, to render the manuscript and the PLOS package:

```bash
Rscript manuscript/make_plos_figures.R
Rscript manuscript/build_plos_submission.R
```

Expect two to four hours on four threads. Almost all of that is
`04_machine_learning.R`: the random hyperparameter search and the ten-fold
out-of-fold predictions run 100+ model fits over 209,000 person-periods. The
search alone takes about 45 minutes and is cached in
`results/ml_tuning.csv`; delete that file, or set `U5M_RETUNE=1`, to force a
fresh search.

### Requirements

R ≥ 4.4 with: `haven`, `dplyr`, `tidyr`, `purrr`, `stringr`, `forcats`,
`data.table`, `sf`, `ggplot2`, `scales`, `viridis`, `patchwork`, `ggrepel`,
`survey`, `broom`, `spdep`, `Matrix`, `INLA`, `xgboost`, `ranger`, `glmnet`,
`yardstick`, `shapviz`, `knitr`, `flextable`, `xml2`. Quarto is needed only
to render the manuscript.

`INLA` is not on CRAN:

```r
install.packages("INLA", repos = c(INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)
```

## Design decisions worth knowing

**Nepal DHS dates are Bikram Sambat.** Every century-month code in these files
— `v008`, `v011`, `b3` — is BS, not Gregorian. Differences between codes (ages,
exposure, birth intervals) are unaffected, but calendar labelling must be
converted or every birth lands about 57 years in the future.
`bs_cmc_to_ad()` in `00_setup.R` does this.

**Person-periods, not a binary outcome.** Each child is expanded into one row
per DHS life-table age segment it entered, with exposure truncated at its
attained age. Fitting a binary "died before five" indicator to a retrospective
birth window instead would credit recently born children with survival through
ages they have not yet reached. Both the life tables and the machine learning
run on the same expanded file, so the descriptive rates, the Bayesian model and
the SHAP decomposition are mutually consistent.

**Ecological belt is not used.** Geography enters the analysis through the
district (the unit of every spatial model), the province and the cluster
altitude. Ecological belt (mountain / hill / Terai) was removed from all
scripts, models, tables and figures.

**DHS code 97 is "not recorded", not "worst category".** Water, sanitation,
cooking fuel and electricity are not collected for a mother who was a visitor
in the sampled household - the same 1,627 children (5.2%) in all four
variables, spread across every wealth quintile including the richest. Coding
them "unimproved"/"polluting" assigned the worst exposure to 5% of the sample
on no evidence; they now carry an explicit `Not recorded` level rather than
being lumped or dropped. Their U5MR (55.9) does sit above that of surveyed
households (46.8), but between the improved (46.2) and unimproved (54.5)
groups, so there is no basis for assigning them either. Because all four
variables are unrecorded for the *same* children, their four indicators are
identical columns: `05` collapses them to one household-survey indicator so the
regression is identifiable, while the tree models keep the three-level coding.

**Ten-year retrospective window.** District-level estimation needs the births.
The cost is recall error, discussed in the manuscript.

**Expected counts are standardised within period.** The scan statistic and the
Bayesian offset use the national age-segment hazard schedule *of the same
period*. Without this the scan simply rediscovers the national decline and
reports the earliest years as a "cluster".

**The tuner picks a large class weight, so probabilities need repair.**
Selecting on precision-recall under a 0.7% event rate drives
`scale_pos_weight` to 20 for both tree models, which multiplies the predicted
odds by 20. `04b` refits Platt scaling on development-set out-of-fold
predictions and applies it to the test set. Discrimination is unchanged (the
transform is monotone); the Brier score for gradient boosting drops from 0.021
to 0.0067 and the calibration slopes move to ~1.

**Validation splits on children, not rows.** All person-periods of a child stay
in the same partition. Leave-one-province-out cross-validation is reported
alongside the random split because neighbouring clusters otherwise leak between
training and test sets and inflate performance.

## Outputs

`results/` holds one CSV per statistic; `results/tables/` holds the display
tables the manuscript renders; `results/key_numbers.csv` holds every figure
quoted in the abstract and results text. The manuscript reads all of these at
render time, so the prose cannot drift away from the analysis.

## Licence

Code: MIT (see `LICENSE`). The NDHS data remain subject to the DHS Program's
terms of use.

## Citation

Acharya J. District-level inequality in under-five mortality in Nepal,
2002–2022: a Bayesian spatio-temporal analysis. Manuscript submitted to
*PLOS ONE*, 2026.
