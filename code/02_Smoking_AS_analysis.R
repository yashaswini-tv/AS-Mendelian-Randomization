# ============================================================
# Smoking Initiation -> Ankylosing Spondylitis
# Two-Sample Mendelian Randomization Analysis
# ============================================================

library(TwoSampleMR)
library(ggplot2)

# Load harmonized smoking initiation-AS dataset
dat <- read.csv("data/Smoking_AS_harmonised.csv")

# Main MR analysis
mr_results <- mr(dat)
print(mr_results)

write.csv(
  mr_results,
  "results/Smoking_AS_MR_results.csv",
  row.names = FALSE
)

# Heterogeneity
heterogeneity <- mr_heterogeneity(dat)
print(heterogeneity)

write.csv(
  heterogeneity,
  "results/Smoking_AS_heterogeneity.csv",
  row.names = FALSE
)

# MR-Egger intercept test
pleiotropy <- mr_pleiotropy_test(dat)
print(pleiotropy)

write.csv(
  pleiotropy,
  "results/Smoking_AS_pleiotropy.csv",
  row.names = FALSE
)

# Single-SNP analysis
single_snp <- mr_singlesnp(dat)

write.csv(
  single_snp,
  "results/Smoking_AS_singleSNP.csv",
  row.names = FALSE
)

# Leave-one-out analysis
loo <- mr_leaveoneout(dat)
print(loo)

write.csv(
  loo,
  "results/Smoking_AS_leaveoneout.csv",
  row.names = FALSE
)

# Instrument strength
dat$F_statistic <- (dat$beta.exposure / dat$se.exposure)^2

summary(dat$F_statistic)
mean(dat$F_statistic, na.rm = TRUE)
min(dat$F_statistic, na.rm = TRUE)
sum(dat$F_statistic < 10, na.rm = TRUE)

# Scatter plot
scatter_plot <- mr_scatter_plot(mr_results, dat)

ggsave(
  "figures/Figure_Smoking_AS_scatter.png",
  plot = scatter_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)

# Leave-one-out plot
loo_plot <- mr_leaveoneout_plot(loo)

ggsave(
  "figures/Figure_Smoking_AS_leave_one_out.png",
  plot = loo_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)

# Funnel plot
funnel_plot <- mr_funnel_plot(single_snp)

ggsave(
  "figures/Figure_Smoking_AS_funnel.png",
  plot = funnel_plot[[1]],
  width = 8,
  height = 6,
  dpi = 300
)
