#!/usr/bin/env Rscript

## Part 2: analyses by individual species
### Which traits explain the number of Wolbachia strains in a given species?
traits_db <- read.table("../European_and_Maghreb_butterfly_trait_data/European_&_Maghreb_Butterfly_Trait_data_v1.2.csv", header = T)
colnames(traits_db)


nuwt <- F

if (nuwt){
  df <- read.csv("../Wolbachia_presence_species_level_infection_conservative.csv", header = T)
  data_type = "infection-conservative"
} else {
  df <- read.csv("../Wolbachia_presence_species_level_NUWT_conservative.csv", header = T)
  data_type = "nuwt-conservative"
}

# Remove Satyrium since we are not working with that genus
df <- df[df$Genus != "Satyrium",]

# Remove species not present in the phylogeny for now (REMEMBER TO ADD THEM BACK LATER!!!)
# missing_spp <- c("ausonia",  "simplonia",  "farinosa",  "rhamni",  "maera",  "panoptes",  "malvoides",  "rosae",  "sertorius",  "lineola",  "polyxena")
# df <- df [!df$Species %in% missing_spp, ]


df$phylo <- paste0(df$Genus, "_", df$Species)

for (speciesnum in 1:length(df$phylo)){
  species <- df$phylo[speciesnum]
  df$voltinism_min[df$phylo == species] <- traits_db$Vol_min[traits_db$Taxon == species]
  df$voltinism_max[df$phylo == species] <- traits_db$Vol_max[traits_db$Taxon == species]
  df$wing_index[df$phylo == species] <- traits_db$WIn[traits_db$Taxon == species]
  df$flight_months_min[df$phylo == species] <- traits_db$FMo_Min[traits_db$Taxon == species]
  df$flight_months_max[df$phylo == species] <- traits_db$FMo_Max[traits_db$Taxon == species]
  df$flight_months_mean[df$phylo == species] <- traits_db$FMo_Average[traits_db$Taxon == species]
}


continuous_vars <- df[, c("voltinism_min", "voltinism_max", "wing_index", "flight_months_min", "flight_months_max", "flight_months_mean")]
print("Pearson correlation")
cor(continuous_vars, method = "pearson")
print("Spearman correlation")
cor(continuous_vars, method = "spearman")

# LIBS
library(ape)

# you may need an ultrametric tree for this? but see: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q3/019095.html
phylo <- read.tree("../../1c_butterfly_phylogeny/phylogeny_full/REROOTED.bfly.phylogeny.full.iqtree.nwk")
str(phylo)


# LIBS
library(phytools)
library(MCMCglmm)
library(stringr)

# MCMCglmm complains about node labels so we save them apart and replace them (they were the bootstrap)
bootstraps <- phylo$node.label
phylo$node.label <- seq(1:length(phylo$node.label))

# change tip labels to match does in the dataframe
phylo$tip.label <- str_to_title(gsub("\\..+", "", phylo$tip.label))

# remove extra numbers present in some species codes
phylo$tip.label <- gsub(pattern = "_[0-9]+", replacement = "", x = phylo$tip.label)

# Compute the inverse branch length matrix
inv.phylo<-inverseA(phylo, nodes='ALL', scale=FALSE)

# Quick check of potential zero inflation, you expect your number of zeroes to be close to the proportion you would expect from a poison with that mean rate
table(df$Number.of.strains==0)
ppois(0, mean(df$Number.of.strains))

hist(df$Number.of.strains)

# in the infection-conservative, there are ~59% of zeroes, so it is somewhat inflated but not much
# in the nuwt-conservative, there are ~45% of zeroes, so it not inflated


# z-transform predictors
df$wing_index_norm <- scale(df$wing_index)
df$voltinism_min_norm <- scale(df$voltinism_min)
df$flight_months_mean_norm <- scale(df$flight_months_mean)

# Specify the prior
# prior <- list(R = list(V = 1, nu = 0))
# The residual variance prior is the one used in the MCMCglmm course notes by Jarrod Hadfield
# The parameter expanded prior for the random effects is recommended in a MCMCglmm tutorial by Pierre de Villemereuil for Poisson models because it is better for variances close to zero, which are common for Poisson because of the exponential term involved in this model (since we have log link function)
prior <- list(R = list(V = 1, nu = 0.002),
              G = list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 1000)
              )
)

# Fit the GLMM model via MCMC
# For both the fixed effects and the random effects, you get the mean of the posterior (which would be your "estimate"), the 95% credible interval (CI), and the effective sample size for that specific term.

# Effective samples sizes should be large (say 1000 or even up to 3000 if you can), as this indicate due to autocorrelation, the number of samples that you have store is as good as if you had that size of independent samples

# For fixed effects you also get a pMCMC value, which is a p-value basically, but may not be so reliable if eff. sample size is small because that means you have a lot of autocorrelated samples (in any case, you are gonna want to increase the eff. sample size anyway!). It is better to evaluate significance with the CI and just looking at the posterior distributions (especially since for the random terms you can't get a p-value)

# For the model to take the phylogenetic effect into account, you MUST include it in the random effects!
# Essentially what is going on is that you model a term of random effect of species which is associated to the phylogeny (is modelled as brownian motion splitting at the nodes, such that two closely related species are more likely to be similar), and also a term of species random effect that is not phhylogenetic, has other sources (so Colias hyale could be more similar to Iphiclides podalirius for example, not to Colias alfacariensis)

# There's also the DIC (deviance information criterion), which is kind of like an AIC but for MCMC apporaches. Smaller DICs should be preferable. However at least in MCMCglmm the DIC value is a poor estimate and it is not recommended to trust it. Since it is based on the same data that was used to build the model it tends to select overfitted models.

# For the random effects/residuals, since the lower 95% CI bound cannot be below zero, a good rule of thumb is that given the effective sample size, it is say 100, then you round the CI to two digits and that should be different from zero.
# If the Poisson is underdispersed just run it as normal and check that both results are consistent
# For the phylogenetic effect, the thing is if it is not ultrametric you cannot scale it, so you getting the mean variance explained by the phylogenetic effect is not straightforward. You can do as follows:
# get average
# If d = distance from root to tip --> get mean_d
# Then mean_phylogenetic_variance_of_the_tips = mean_d * variance estimate from model (mean of the posterior of the random effect for phylo)
# Then proportion of variance that is phylogenetic = mean_phylogenetic_variance_of_the_tips/mean_phylogenetic_variance_of_the_tips+residual variance (to get the total)
model_pois <- MCMCglmm(Number.of.strains ~ wing_index_norm * voltinism_min_norm * flight_months_mean_norm,
                       data = df,
                       random = ~ phylo,
                       ginverse = list(phylo=inv.phylo$Ainv),
                       family = "poisson",
                       prior = prior,
                       thin = 450,
                       burnin = 3000,
                       nitt = 4503000, pl = T)

simulation <- simulate.MCMCglmm(model_pois, nsim = 1000)

model_gaus <- MCMCglmm(Number.of.strains ~ wing_index_norm * voltinism_min_norm * flight_months_mean_norm,
                       data = df,
                       random = ~ phylo,
                       ginverse = list(phylo=inv.phylo$Ainv),
                       family = "gaussian",
                       prior = prior,
                       thin = 150,
                       burnin = 3000,
                       nitt = 1503000, pl = T)

# checking how sensible the number of outliers in the model are
strain_num_vect <- c()
cumu_prob_vect <- c()

for (i in 0:100){
  strain_num_vect[length(strain_num_vect)+1] <- i
  cumu_prob_vect[length(cumu_prob_vect)+1] <- sum(simulation >= i)/sum(simulation >= 0)
}

plot(strain_num_vect, cumu_prob_vect)
# the probability of having fairly large number of strains is not that low, so this is a caveat of fitting the Poisson in this case, but we have no clear upper bound either so it is probably OK to fit it, it is still rare outliers, but important to keep in mind.

# write to files
  unlink(paste0("mcmcglmm_part2_model_", data_type, ".txt"))
  sink(file = paste0("mcmcglmm_part2_model_", data_type, ".txt"))
  print("Poisson model")
  print(summary(model_pois))
  print("Gaussian model")
  print(summary(model_gaus))
  sink()
  
  pdf(paste0("mcmcglmm_part2_model_", data_type, "_fixed_effects.pdf"),
      width = 7, height = 14)
  par(mfrow = c(8,2), mar = c(2, 2, 2, 2))
  for(n in 1:8){
    # svg(file = paste0("mcmcglmm_part2_model_infection-conservative_fixed_effects_", as.character(n), ".svg"),
        # width = 12, height = 4)
    plot(seq(from = 3001, to = 4503000, by = 450), model_pois$Sol[,n], type = "l",
         main = paste0("Trace of ", colnames(model_pois$Sol)[n]),
         xlab = "Iterations",
         ylab = NA,
         cex.main = 0.8)
    # plot(model_pois$Sol[, n], cex.main = 0.8)
    plot(density(model_pois$Sol[,n]),
         main = paste0("Density of ", colnames(model_pois$Sol)[n]),
         xlab = "N = 10,000",
         cex.main = 0.8)
  }
  dev.off()
  pdf(file = paste0("mcmcglmm_part2_model_", data_type, "_random_effects.pdf"))
  plot(model_pois$VCV)
  dev.off()
  pdf(file = paste0("mcmcglmm_part2_model_", data_type, "_histogram.pdf"))
  hist(colSums(simulation==0))
  abline(v = length(df$Number.of.strains[df$Number.of.strains==0]), col = "red")
  dev.off()
  pdf(paste0("mcmcglmm_part2_model_", data_type, "_fixed_effects_gaussian.pdf"))
  plot(model_gaus$Sol, cex.main = 0.8)
  plot(model_gaus$VCV)
  dev.off()
  
  # save.image(file = paste0("mcmcglmm_part2_model_", data_type, ".RData"))

  # unlink("mcmcglmm_part2_model_nuwt-conservative.txt")
  # sink(file = "mcmcglmm_part2_model_nuwt-conservative.txt")
  # print("Poisson model")
  # print(summary(model_pois))
  # print("Gaussian model")
  # print(summary(model_gaus))
  # sink()
  # 
  # svg(file = "mcmcglmm_part2_model_nuwt-conservative_fixed_effects_main.svg", 
  #     width = 12, height = 4)
  # plot(model_pois$Sol[, c(1:4)], cex.main = 0.8)
  # dev.off()
  # svg(file = "mcmcglmm_part2_model_nuwt-conservative_fixed_effects_int.svg")
  # plot(model_pois$Sol[, c(5:8)], cex.main = 0.8)
  # dev.off()
  # svg(file = "mcmcglmm_part2_model_nuwt-conservative_random_effects.svg")
  # plot(model_pois$VCV)
  # dev.off()
  # svg(file = "mcmcglmm_part2_model_nuwt-conservative_histogram.svg")
  # hist(colSums(simulation==0))
  # abline(v = length(df$Number.of.strains[df$Number.of.strains==0]), col = "red")
  # dev.off()
  # plot(model_gaus$Sol)
  # plot(model_gaus$VCV)
  # dev.off()
  # 
  # # save.image(file = "mcmcglmm_part2_model_nuwt-conservative.RData")
