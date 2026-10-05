# Reproducibility

## Software

The Mendelian randomization analyses were conducted in R using the `TwoSampleMR` package. MR-PRESSO was used as an additional sensitivity analysis for the BMI analysis.

## Analysis

Separate analysis scripts are provided for:

- BMI → ankylosing spondylitis
- Smoking initiation → ankylosing spondylitis

The scripts are available in the `code/` directory.

## Data Availability

This project uses summary-level genome-wide association study (GWAS) data.

Large raw GWAS summary-statistic files are not included in this repository because of their size. The repository contains analysis code, selected result files, and figures required to document the analyses.

## Reproducing the Analysis

The analysis scripts expect harmonized exposure-outcome datasets before the MR analyses are run. Local file paths may need to be adjusted according to the user's computing environment.

## Important Note

The repository is intended to document the statistical analysis and improve transparency and reproducibility. Raw GWAS datasets should be obtained from their original data providers.
