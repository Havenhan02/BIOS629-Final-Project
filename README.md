# BIOS629 Final Project: Clinical vs. Genomic Predictors of Overall Survival in NSCLC

**Author:** Wenjie Han  
**Course:** BIOS 629 - Winter 2026  
**University of Michigan, Department of Biostatistics**

## Overview

This project investigates whether genomic features (Fraction Genome Altered, Mutation Count, TMB) can improve prediction of overall survival in Non-Small Cell Lung Cancer (NSCLC) beyond standard clinical variables (stage, age, sex, smoking history).

## Data Source

- **Dataset:** MSK-CHORD (Memorial Sloan Kettering Cancer Center, Nature 2024)
- **Population:** 7,809 NSCLC patients
- **Source:** [cBioPortal - MSK-CHORD 2024](https://www.cbioportal.org/study/summary?id=msk_chord_2024)

## Reproducibility

The main analysis file `final_project.Rmd` can be knit directly to reproduce all results. Data is read from this repository via GitHub raw URL.

### Requirements

- R (≥ 4.0)
- R packages: `tidyverse`, `survival`, `survminer`, `randomForestSRC`, `survcomp`, `caret`, `gridExtra`

### Running the Analysis

1. Clone this repository
2. Open `final_project.Rmd` in RStudio
3. Click **Knit** to generate the full HTML report

All numerical results from deterministic models (Cox PH) will be exactly reproduced. Results from random forest (RSF) may vary slightly across R versions but should be similar with `set.seed(629)`.

## Repository Structure

```
├── final_project.Rmd          # Main analysis (R Markdown)
├── run_analysis.R             # R script version
├── gen_risk.R                 # Risk group generation script
├── generate_km_plots.R        # KM plot generation script
├── generate_plots.R           # Other plot generation scripts
├── dataset/
│   └── msk_chord_2024_NSCLC_data.csv   # NSCLC subset data
└── final_output/
    ├── fig1_cv_boxplot.png              # Cross-validation C-index boxplot
    ├── fig2_cv_barplot.png              # CV performance bar plot
    ├── fig3_variable_importance.png     # Variable importance plot
    ├── fig4_km_overall.png              # Overall KM survival curve
    ├── fig5_km_stage.png                # KM by cancer stage
    ├── fig6_km_risk_groups.png          # KM by risk groups
    ├── risk_group_summary.csv           # Risk group statistics
    ├── analysis_risk.rds                # RSF risk analysis object
    ├── cv_results.rds                   # Cross-validation results
    ├── cv_summary.rds                   # CV summary statistics
    └── importance_df.rds                # Variable importance data
```

## Key Results

| Model | 5-Fold CV C-index | Full-Data C-index |
|-------|-------------------|-------------------|
| Clinical Cox PH | 0.668 ± 0.015 | 0.669 |
| Clinical + Genomic Cox PH | 0.678 ± 0.015 | 0.679 |
| Random Survival Forest | 0.681 ± 0.014 | 0.691 |

**Top Predictors:** Stage > Fraction Genome Altered > Mutation Count > TMB > Age > Smoking > Sex

**Risk Stratification (Tertiles):**
- Low Risk: 34.3 mo median OS, 28% event rate
- Intermediate Risk: 22.2 mo median OS, 53.2% event rate  
- High Risk: 13.5 mo median OS, 73.1% event rate
