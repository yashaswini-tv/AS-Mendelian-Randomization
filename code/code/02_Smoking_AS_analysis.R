# ============================================================
# Smoking Initiation -> Ankylosing Spondylitis
# Two-Sample Mendelian Randomization Analysis
# ============================================================

library(TwoSampleMR)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load harmonized smoking initiation-AS data
# ------------------------------------------------------------

# This file contains the 19 SNPs retained after harmonization
dat <- read.csv("data/Smoking_AS_harmonised.csv")

dim(dat)
head(dat)

# ------------------------------------------------------------
# 2. Mendelian Randomization Analysis
# ------------------------------------------------------------

mr_results <- mr(dat)

print(mr_results)

write.csv(
  mr_results,
  "results/Smoking_AS_MR_results.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 3. Heterogeneity
# ------------------------------------------------------------

heterogeneity <- mr_heterogeneity(dat)

print(heterogeneity)

write.csv(
  heterogeneity,
  "results/Smoking_AS_heterogeneity.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 4. Horizontal Pleiotropy
# ------------------------------------------------------------

pleiotropy <- mr_pleiotropy_test(dat)

print(pleiotropy)

write.csv(
  pleiotropy,
  "results/Smoking_AS_pleiotropy.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Single-SNP Analysis
# ------------------------------------------------------------

single_snp <- mr_singlesnp(dat)

write.csv(
  single_snp,
  "results/Smoking_AS_singleSNP.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Leave-One-Out Analysis
# ------------------------------------------------------------

loo <- mr_leaveoneout(dat)

print(loo)

write.csv(
  loo,
  "results/Smoking_AS_leaveoneout.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Instrument Strength
# ------------------------------------------------------------

dat$F_statistic <- (dat$beta.exposure / dat$se.exposure)^2

summary(dat$F_statistic)

mean(dat$F_statistic, na.rm = TRUE)
min(dat$F_statistic, na.rm = TRUE)
sum(dat$F_statistic < 10, na.rm = TRUE)

# ------------------------------------------------------------
# 8. MR Scatter Plot
# ------------------------------------------------------------

scatter_plot <- mr_scatter_plot(mr_results, dat)

ggsave(
  "figures/Figure2_Smoking_scatter.png",
  plot = scatter_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 9. Leave-One-Out Plot
# ------------------------------------------------------------

loo_plot <- mr_leaveoneout_plot(loo)

ggsave(
  "figures/Figure3_Smoking_leave_one_out.png",
  plot = loo_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Funnel Plot
# ------------------------------------------------------------

funnel_plot <- mr_funnel_plot(single_snp)

ggsave(
  "figures/Figure4_Smoking_funnel.png",
  plot = funnel_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# End of smoking initiation-AS analysis
# ------------------------------------------------------------
