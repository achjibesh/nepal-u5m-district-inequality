---
output:
  word_document: default
  html_document: default
---
# S2 Checklist. TRIPOD+AI checklist

**Manuscript:** District-level inequality in under-five mortality in Nepal, 2002–2022: a Bayesian spatio-temporal analysis

TRIPOD+AI statement (Collins et al., *BMJ* 2024;385:e078378), applied to the
prediction-model component: elastic-net logistic regression, random forest,
gradient boosting and a stacked ensemble for death within a DHS age segment.
Locations refer to manuscript sections; page and line numbers should be added
from the final typeset file.

| Section | Item | Location in manuscript |
|---|---|---|
| **Title** | 1. Identify the study as developing or evaluating a prediction model, the target population and outcome | Abstract, Methods and Results (the title names the primary Bayesian analysis; the prediction-model component, population and outcome are identified in the abstract) |
| **Abstract** | 2. Structured summary | Abstract |
| **Introduction** | 3a. Healthcare context and rationale, including existing models | Introduction, paragraph 4 |
| | 3b. Target population and intended use | Introduction, final paragraph; statement of intended use below |
| | 3c. Known inequalities between groups | Introduction, paragraph 2 (caste/ethnic and provincial inequalities) |
| | 4. Objectives, including development and/or validation | Introduction, final paragraph (aim iv) |
| **Methods** | 5a. Data sources, separately for development and evaluation | Materials and methods: Data; Machine learning (70/30 child split) |
| | 5b. Dates of data collection and follow-up | Materials and methods: Data; Study population and outcome |
| | 6a. Setting | Materials and methods: Data; Linking clusters to districts |
| | 6b. Eligibility criteria | Materials and methods: Study population and outcome; Table S1 in S1 File |
| | 6c. Treatments received, if relevant | Not applicable |
| | 7. Data preparation and quality checks | Materials and methods: Study population and outcome; Covariates; Section S1 in S1 File |
| | 8a. Outcome definition and how and when it was assessed | Materials and methods: Study population and outcome |
| | 8b. Blinding of outcome assessment | Not applicable (outcome from maternal birth-history report) |
| | 8c. Outcome variability across groups | Table 1; Table 2 |
| | 9a. Predictors, including rationale for choice | Materials and methods: Covariates (Mosley–Chen framework) |
| | 9b. Predictor definition and timing of measurement | Materials and methods: Covariates; Discussion: Strengths and limitations (measurement at interview) |
| | 9c. Predictors excluded and why | Materials and methods: Covariates (children ever born and household size excluded as outcome-dependent; maternal care variables unavailable for most of the window) |
| | 10. Sample size and number of events | Results: Sample and national trends; Table S1 in S1 File |
| | 11. Missing data | Materials and methods: Covariates (first births; DHS code 97 retained as a category) |
| | 12a. Data use: development and evaluation | Materials and methods: Machine learning |
| | 12b. Predictor handling | Materials and methods: Machine learning (one encoding per variable, continuous where available) |
| | 12c. Model type, building steps and internal validation | Materials and methods: Machine learning (random search scored on five of ten child-grouped cross-validation folds, ten-fold out-of-fold stacking, single-use test set); Tables S9 and S10 in S1 File |
| | 12d. Class imbalance | Materials and methods: Machine learning (class weights tuned; AUPRC used; Platt recalibration); Table S12 in S1 File |
| | 12e. Fairness and group comparisons | Materials and methods: Machine learning (leave-one-province-out validation; province-stratified SHAP); caste/ethnicity included as a predictor |
| | 12f. Model output | Materials and methods: Machine learning (segment hazard; child-level probability of death before age five) |
| | 12g. Heterogeneity across clusters/settings | Materials and methods: Machine learning (leave-one-province-out validation) |
| | 13. Performance measures and rationale | Materials and methods: Machine learning (AUROC, AUPRC, Brier score, calibration intercept and slope, decile calibration) |
| | 14. Model updating or recalibration | Materials and methods: Machine learning (Platt scaling of tree models) |
| | 15. Interpretability methods | Materials and methods: Machine learning (TreeSHAP; province-stratified SHAP; partial dependence) |
| | 16. Hardware and software | Materials and methods: Software and ethics; Table S20 in S1 File |
| **Open science** | 17a. Funding | Funding |
| | 17b. Conflicts of interest | Declaration of competing interest |
| | 17c. Protocol | Not prepared (secondary analysis) |
| | 17d. Registration | Not registered (secondary analysis of public-use data) |
| | 18a. Data availability | Data availability |
| | 18b. Code availability | Data availability |
| **Patient and public involvement** | 19. Details of involvement | Not applicable (secondary analysis of anonymised survey data) |
| **Results** | 20a. Participant flow | Table S1 and Fig S2 in S1 File |
| | 20b. Characteristics of participants | Table 1; Results: Sample and national trends |
| | 20c. Comparison of development and evaluation data | Materials and methods: Machine learning (split stratified by outcome and survey round) |
| | 21. Number of participants and events in each analysis | Results: Sample and national trends; Table S11 in S1 File |
| | 22. Full model specification | Table 6 (regression counterpart); fitted model objects available with the code |
| | 23a. Performance estimates with confidence intervals | Table 7; Fig 5 |
| | 23b. Performance across groups | Fig 5D; Table S11 in S1 File |
| | 24. Model updating results | Tables S12 and S13 in S1 File; Fig 5C |
| | 25. Interpretability results | Figs 6 and 7; Table S14 and Fig S7 in S1 File |
| **Discussion** | 26. Interpretation, including comparison with existing models | Discussion: Machine learning |
| | 27a. Handling of poor-quality or unavailable input data in use | Not applicable (model not intended for deployment) |
| | 27b. User interaction and expertise required | Not applicable (model not intended for deployment) |
| | 27c. Future research and generalisability | Discussion: Residual district variation; Strengths and limitations |
| | 28. Limitations, including generalisability | Discussion: Strengths and limitations |

## Statement of intended use

The prediction models are not intended for individual clinical decisions or for
screening individual children. They are used to assess how well measured
characteristics predict child death, whether machine learning adds to penalised
regression, whether a model transfers across provinces, and whether the
contribution of predictors differs geographically. Discrimination and
calibration are reported so that these interpretive conclusions can be judged
against the models' actual predictive performance.
