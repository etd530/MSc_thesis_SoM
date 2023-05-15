#!/usr/bin/env Rscript

# Set boolean variable to run either one analysis or the other
nuwt = T

if (nuwt){
  df <- read.csv("../Species_pairs_traits_infection_conservative.csv", header = T)
} else {
  df <- read.csv("../Species_pairs_traits_NUWT_conservative.csv", header = T)
}

df$Genus <- gsub("*", "", df$Genus, fixed = T)

df <- df[df$Genus != "Fabriciana" & df$Genus != "Erebia",]

df$Split_time_mean <- as.double(gsub("\\(.+\\)", "", df$Split.time..MYA.))

library(plyr)
df$Gen.y.1.sp.1 <- as.integer(laply(.data = strsplit(df$Gen.y.1.sp.1, split = "–", fixed = T), .fun = function(x) x[[length(x)]]))
df$Gen.y.1.sp.2 <- as.integer(laply(.data = strsplit(df$Gen.y.1.sp.2, split = "–", fixed = T), .fun = function(x) x[[length(x)]]))

for (rownum in 1:nrow(df)){
  df$gen.y.1.mean[rownum] <- mean(c(df$Gen.y.1.sp.1[rownum], df$Gen.y.1.sp.2[rownum]))
}

# The residual prior is important if we have random effects --> revise MCMCglmm documentation if you add them, but otherwise I think it is fine?
prior <- list(R = list(V = 1, fix = 1))
# B = list(mu = c(0, 0), V = diag(2) * + (1 + pi^2/3)),

# build subset data
df_subset <- df[!(is.na(df$Presence.of.shared.strains.in.the.pair)),]
df_subset$Presence.of.shared.strains.in.the.pair[df_subset$Presence.of.shared.strains.in.the.pair==1] <- "Y"
df_subset$Presence.of.shared.strains.in.the.pair[df_subset$Presence.of.shared.strains.in.the.pair==0] <- "N"
df_subset$Presence.of.shared.strains.in.the.pair <- as.factor(df_subset$Presence.of.shared.strains.in.the.pair)

# Z-transform predictors
df_subset$Degree.of.sympatry_norm <- scale(df_subset$Degree.of.sympatry)
df_subset$Split_time_mean_norm <- scale(df_subset$Split_time_mean)

library(MCMCglmm)

# Fit the GLMM model via MCMC
# you want to fit categorical, multinomial2 wants the count of each type of outcome, not just the true/false values of the outcome, that's why it complains, it needs to response variables
# Also we fit a full model since there may be interactions and you don't want biased estimates
model <- MCMCglmm(Presence.of.shared.strains.in.the.pair ~ Split_time_mean_norm * Degree.of.sympatry_norm -1,
                  data = df_subset,
                  family = "threshold",
                  thin = 500,
                  burnin = 3000,
                  nitt = 5003000,
                  prior = prior)

# write to files
if (nuwt){
  unlink("mcmcglmm_logistic_model_results_infection-conservative_normalised_split_time_mya.txt")
  sink(file = "mcmcglmm_logistic_model_results_infection-conservative_normalised_split_time_mya.txt")
  print(summary(model))
  sink()
  
  pdf(file = "mcmcglmm_logistic_model_results_infection-conservative_normalised_split_time_mya.pdf")
  plot(model$Sol)
  plot(model$VCV)
  dev.off()
  
  save.image(file = "mcmcglmm_logistic_model_results_infection_conservative_normalised_split_time_mya.RData")
  
  print("results saved")
} else {
  unlink("mcmcglmm_logistic_model_results_nuwt-conservative_normalised_split_time_mya.txt")
  sink(file = "mcmcglmm_logistic_model_results_nuwt-conservative_normalised_split_time_mya.txt")
  print(summary(model))
  sink()
  
  pdf(file = "mcmcglmm_logistic_model_results_nuwt-conservative_normalised_split_time_mya.pdf")
  plot(model$Sol)
  plot(model$VCV)
  dev.off()
  
  save.image(file = "mcmcglmm_logistic_model_results_nuwt-conservative_normalised_split_time_mya.RData")
}
