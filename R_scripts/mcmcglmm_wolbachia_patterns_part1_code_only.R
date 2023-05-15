#!/usr/bin/env Rscript

# Set boolean variable to run either one analysis or the other
nuwt = F

if (nuwt){
  df <- read.csv("../Species_pairs_traits_infection_conservative.csv", header = T)
  data_type = "infection-conservative"
} else {
  df <- read.csv("../Species_pairs_traits_NUWT_conservative.csv", header = T)
  data_type = "nuwt-conservative"
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
df_subset$split.time..gen.norm <- scale(df_subset$split.time..gen.)

library(MCMCglmm)

# Fit the GLMM model via MCMC
# you want to fit categorical, multinomial2 wants the count of each type of outcome, not just the true/false values of the outcome, that's why it complains, it needs to response variables
# Also we fit a full model since there may be interactions and you don't want biased estimates
model <- MCMCglmm(Presence.of.shared.strains.in.the.pair ~ split.time..gen.norm * Degree.of.sympatry_norm -1,
                  data = df_subset,
                  family = "threshold",
                  thin = 50,
                  burnin = 3000,
                  nitt = 503000,
                  prior = prior)

# write to files
  unlink(paste0("mcmcglmm_logistic_model_results_", data_type, "_normalised.txt"))
  sink(file = paste0("mcmcglmm_logistic_model_results_", data_type, "_normalised.txt"))
  print(summary(model))
  sink()
  
  pdf(file = paste0("mcmcglmm_logistic_model_results_", data_type, "_normalised.pdf"))
  par(mfrow = c(3,2), mar = c(2, 2, 2, 2))
  for(n in 1:3){
    plot(seq(from = 3001, to = 503000, by = 50), model$Sol[,n], type = "l",
         main = paste0("Trace of ", colnames(model$Sol)[n]),
         xlab = "Iterations",
         ylab = NA,
         cex.main = 0.8)

    plot(density(model$Sol[,n]),
         main = paste0("Density of ", colnames(model$Sol)[n]),
         xlab = "N = 10,000",
         cex.main = 0.8)
  }
  # plot(model$Sol)
  # plot(model$VCV)
  dev.off()
  
  # save.image(file = "mcmcglmm_logistic_model_results_infection_conservative_normalised.RData")
  
  # print("results saved")
  # 
  # unlink("mcmcglmm_logistic_model_results_nuwt-conservative_normalised.txt")
  # sink(file = "mcmcglmm_logistic_model_results_nuwt-conservative_normalised.txt")
  # print(summary(model))
  # sink()
  # 
  # pdf(file = "mcmcglmm_logistic_model_results_nuwt-conservative_normalised.pdf")
  # plot(model$Sol)
  # plot(model$VCV)
  # dev.off()
  # 
  # save.image(file = "mcmcglmm_logistic_model_results_nuwt-conservative_normalised.RData")
