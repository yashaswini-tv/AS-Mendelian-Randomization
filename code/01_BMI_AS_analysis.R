# ============================================================
# BMI -> Ankylosing Spondylitis
# Two-Sample Mendelian Randomization Analysis
# ============================================================

# Packages
library(TwoSampleMR)
library(MRPRESSO)

# ------------------------------------------------------------
# 1. Load harmonized BMI-AS data
# ------------------------------------------------------------

# Change this path to the location of BMI_AS_harmonised.csv
dat <- read.csv("data/BMI_AS_harmonised.csv")

# Check data
dim(dat)
head(dat)

# ------------------------------------------------------------
# 2. Mendelian Randomization Analysis
# ------------------------------------------------------------

mr_results <- mr(dat)

print(mr_results)

write.csv(
  mr_results,
  "results/BMI_AS_MR_results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Heterogeneity Analysis
# ------------------------------------------------------------

heterogeneity <- mr_heterogeneity(dat)

print(heterogeneity)

write.csv(
  heterogeneity,
  "results/BMI_AS_heterogeneity.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Horizontal Pleiotropy
# ------------------------------------------------------------

pleiotropy <- mr_pleiotropy_test(dat)

print(pleiotropy)

write.csv(
  pleiotropy,
  "results/BMI_AS_pleiotropy.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Single-SNP Analysis
# ------------------------------------------------------------

single_snp <- mr_singlesnp(dat)

write.csv(
  single_snp,
  "results/BMI_AS_singleSNP.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Leave-One-Out Analysis
# ------------------------------------------------------------

loo <- mr_leaveoneout(dat)

write.csv(
  loo,
  "results/BMI_AS_leaveoneout.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Funnel Plot
# ------------------------------------------------------------

funnel_plot <- mr_funnel_plot(single_snp)

ggsave(
  "figures/Figure1_BMI_funnel.png",
  plot = funnel_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Instrument Strength
# ------------------------------------------------------------

dat$F_statistic <- (dat$beta.exposure / dat$se.exposure)^2

summary(dat$F_statistic)

mean(dat$F_statistic, na.rm = TRUE)
min(dat$F_statistic, na.rm = TRUE)
sum(dat$F_statistic < 10, na.rm = TRUE)

# ------------------------------------------------------------
# 9. MR-PRESSO
# ------------------------------------------------------------

set.seed(123)

presso_results <- mr_presso(
  BetaOutcome = "beta.outcome",
  BetaExposure = "beta.exposure",
  SdOutcome = "se.outcome",
  SdExposure = "se.exposure",
  OUTLIERtest = TRUE,
  DISTORTIONtest = TRUE,
  data = dat,
  NbDistribution = 1000,
  SignifThreshold = 0.05
)

print(presso_results)

saveRDS(
  presso_results,
  "results/BMI_AS_MRPRESSO_results.rds"
)

# ------------------------------------------------------------
# 10. Sensitivity Analysis:
#     Exclude Two Extreme BMI Variants
# ------------------------------------------------------------

extreme_snps <- c("rs4713436", "rs4151664")

dat_sensitivity <- dat[
  !(dat$SNP %in% extreme_snps),
]

# Repeat MR analysis
sensitivity_mr <- mr(dat_sensitivity)

print(sensitivity_mr)

write.csv(
  sensitivity_mr,
  "results/BMI_AS_outlier_sensitivity_MR.csv",
  row.names = FALSE
)

# Repeat heterogeneity analysis
sensitivity_heterogeneity <- mr_heterogeneity(dat_sensitivity)

print(sensitivity_heterogeneity)

write.csv(
  sensitivity_heterogeneity,
  "results/BMI_AS_outlier_sensitivity_heterogeneity.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# End of BMI-AS analysis
# ------------------------------------------------------------
