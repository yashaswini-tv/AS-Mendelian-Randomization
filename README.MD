# Causal Effects of Body Mass Index and Smoking Initiation on Ankylosing Spondylitis

## Overview

This repository contains the analysis code and results for a two-sample Mendelian randomization (MR) study investigating the potential causal effects of body mass index (BMI) and smoking initiation on ankylosing spondylitis (AS).

## Methods

The following Mendelian randomization methods were used:

- Inverse-variance weighted (IVW)
- MR-Egger
- Weighted median
- Simple mode
- Weighted mode

Sensitivity analyses included Cochran's Q heterogeneity test, MR-Egger intercept, leave-one-out analysis, MR-PRESSO, and BMI outlier sensitivity analysis.

## Main Results

No statistically significant evidence supporting causal effects of genetically predicted BMI or smoking initiation on ankylosing spondylitis was identified.

### BMI
- SNPs: 491
- IVW beta: -0.095
- SE: 0.113
- P-value: 0.401

### Smoking Initiation
- SNPs: 19
- IVW beta: -0.276
- SE: 0.260
- P-value: 0.289

## Software

Analyses were conducted in R using the TwoSampleMR package. MR-PRESSO was additionally used for sensitivity analysis.

## Data

This study uses summary-level genome-wide association study (GWAS) data. Large raw GWAS summary-statistic files are not included in this repository.
