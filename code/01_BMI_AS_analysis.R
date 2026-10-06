# ============================================================
# BMI -> Ankylosing Spondylitis
# Two-Sample Mendelian Randomization Analysis
# ============================================================

library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)

# Load harmonized BMI-AS dataset
dat <- read.csv("data/BMI_AS_harmonised.csv")

# Main MR analysis
mr_results <- mr(dat)
print(mr_results)

write.csv(
  mr_results,
  "results/BMI_AS_MR_results.csv",
  row.names = FALSE
)

# Heterogeneity
heterogeneity <- mr_heterogeneity(dat)
print(heterogeneity)

write.csv(
  heterogeneity,
  "results/BMI_AS_heterogeneity.csv",
  row.names = FALSE
)

# MR-Egger intercept test
pleiotropy <- mr_pleiotropy_test(dat)
print(pleiotropy)

write.csv(
  pleiotropy,
  "results/BMI_AS_pleiotropy.csv",
  row.names = FALSE
)

# Single-SNP analysis
single_snp <- mr_singlesnp(dat)

write.csv(
  single_snp,
  "results/BMI_AS_singleSNP.csv",
  row.names = FALSE
)

# Leave-one-out analysis
loo <- mr_leaveoneout(dat)

write.csv(
  loo,
  "results/BMI_AS_leaveoneout.csv",
  row.names = FALSE
)

# Funnel plot
funnel_plot <- mr_funnel_plot(single_snp)

ggsave(
  "figures/Figure_5_BMI_AS_funnel.png",
  plot = funnel_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)

# Instrument strength
dat$F_statistic <- (dat$beta.exposure / dat$se.exposure)^2

summary(dat$F_statistic)
mean(dat$F_statistic, na.rm = TRUE)
min(dat$F_statistic, na.rm = TRUE)
sum(dat$F_statistic < 10, na.rm = TRUE)

# MR-PRESSO
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

# Sensitivity analysis excluding two extreme variants
extreme_snps <- c("rs4713436", "rs4151664")

dat_sensitivity <- dat[
  !(dat$SNP %in% extreme_snps),
]

sensitivity_mr <- mr(dat_sensitivity)
print(sensitivity_mr)

write.csv(
  sensitivity_mr,
  "results/BMI_AS_outlier_sensitivity_MR.csv",
  row.names = FALSE
)

sensitivity_heterogeneity <- mr_heterogeneity(dat_sensitivity)
print(sensitivity_heterogeneity)

write.csv(
  sensitivity_heterogeneity,
  "results/BMI_AS_outlier_sensitivity_heterogeneity.csv",
  row.names = FALSE
)
