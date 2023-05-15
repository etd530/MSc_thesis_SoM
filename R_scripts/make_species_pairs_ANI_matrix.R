#!/usr/bin/env Rscript

# Script to test if mean ANI between Wolbachia strains present in sister pairs is
# higher or lower than when comparing random species pairs.

# The rationale is that if Wolbachia tends to codiverge with its host early in the
# speciation, then you would expect ANI between strains present in recently-diverged
# sister taxa to be lower than when comparing strains present in other taxa, as there
# is a trend for these strains to share a recent last common ancestor


#### LIBRARIES ####
library(tidyverse)

#### RAW DATA ####

# Since we have two versions of the data we run the script twice with different dataframes

nuwt <- F

if(nuwt) {
  df <- read.csv("../Samples_all_wolbachia_presence_absence_infection_conservative.csv")
} else {
  df <- read.csv("../Samples_all_wolbachia_presence_absence_NUWT_conservative.csv")
}

ANI_matrix <- read.table("/home/etd530/Documents/TFM_holobiome/13_reseq_datasets_mapping/genomes_all_good_quality/wg_ani/fastani.all_vs_all.txt")
colnames(ANI_matrix) <- c("Strain1", "Strain2", "ANI", "Fragments_aligned", "Fragments_total")

#### MAIN ####
df %>% unite(Genus, Species, col = "Species_binomial", sep = " ", remove = F) -> df

# Discard species that are not sister taxa with any of the ones we have
df <- df[!(df$Species_binomial %in% c("Brenthis hecate",
                                      "Euchloe simplonia",
                                      "Gonepteryx farinosa",
                                      "Lasiommata maera",
                                      "Melanargia ines",
                                      "Pseudophilotes vicrama",
                                      "Spialia rosae",
                                      "Thymelicus lineola")),]

# Discard uninfected samples
df <- df[!is.na(df$Closest.reference.Wolbachia),]
tail(df)

# Keep only genera where both species of the pair are infected
infected_genera <- c()
for (genus in unique(df$Genus)) {
  species_num <- length(unique(df$Species_binomial[df$Genus == genus]))
  if(species_num>1) {infected_genera[length(infected_genera)+1] <- genus}
}
df <- df[df$Genus %in% infected_genera,]

# Since coinfections are listed in a single cell separated by commas,
# we duplicate entries of samples coinfected and put one strain in each entry
# to ease the counting of strains
for (rownum in 1:nrow(df)){
  strains <- df$Closest.reference.Wolbachia[rownum]
  if (grepl(",", strains, fixed = T)){
    strains_splitted <- strsplit(strains, ", ")[[1]]
    for (strain_num in 1:length(strains_splitted)){
      if(strain_num == 1){
        df[rownum, "Closest.reference.Wolbachia"] <- strains_splitted[strain_num]
      } else {
        df[nrow(df) +1, ] <- df[rownum, ]
        df[nrow(df), "Closest.reference.Wolbachia"] <- strains_splitted[strain_num]
      }
    }
  }
}
tail(df)

# Erebia specimens mapped equally well to euryale_EE_932 and ligea_RO_EL_949, they always mapped slightly better to
# euryale so we consider them that strain for this analysis
df$Closest.reference.Wolbachia <- gsub("w.erebia_euryale.EE_932.v1.fna/w.erebia_ligea.RO_EL_949.v1.fna", 
                                       "w.erebia_euryale.EE_932.v1.fna", df$Closest.reference.Wolbachia, fixed = T)

# For P. eros we also have two strains with a slash, we will consider them the P. eros strain since it gets substantially higher coverages
df$Closest.reference.Wolbachia <- gsub("w.polyommatus_eros.PE_1417.v2.fna/w.colias_crocea.Lep_ilColCroc_2.fna",
                                       "w.polyommatus_eros.PE_1417.v2.fna", 
                                       df$Closest.reference.Wolbachia, fixed = T)

head(df)

# Build matrix of N x N
species_list = unique(df$Species_binomial)
species_num = length(species_list)
pairs_matrix <- as.data.frame(matrix(nrow = species_num, ncol = species_num))
colnames(pairs_matrix) = species_list
rownames(pairs_matrix) = species_list


for (species1 in species_list){
  subset1 <- df[df$Species_binomial == species1,]
  for (species2 in species_list){
    ANI_list <- c()
    subset2 <- df[df$Species_binomial == species2,]
    for (samplenum1 in 1:nrow(subset1)){
      strain1 <- subset1$Closest.reference.Wolbachia[samplenum1]
      for (samplenum2 in 1:nrow(subset2)){
        strain2 <- subset2$Closest.reference.Wolbachia[samplenum2]
        ANI_list[length(ANI_list)+1] <- ANI_matrix$ANI[ANI_matrix$Strain1 == strain1 & ANI_matrix$Strain2 == strain2]
        ANI_list[length(ANI_list)+1] <- ANI_matrix$ANI[ANI_matrix$Strain1 == strain2 & ANI_matrix$Strain2 == strain1]
      }
    }
    pairs_matrix[species1, species2] <- mean(ANI_list)
  }
}

# Build fake sets of sister pairs and calculate average number of strains shared in each pair
# Vector to store mean number of strains shared in each fake dataset
mean_ani_list <- c()

for (replicate in 1:19999){
  # Create vector to store ANI in each of the fake pairs
  ani_list <- c()
  # get full species list
  species_list = unique(df$Species_binomial)
  while (length(species_list) > 0){
    # Get first species and drop it from the species list
    num <- round(runif(1, min = 1, max = length(species_list)))
    sister1 <- species_list[num]
    species_list <- species_list[-num]
    
    # Get second species and drop it from the species list  
    num <- round(runif(1, min = 1, max = length(species_list)))
    sister2 <- species_list[num]
    species_list <- species_list[-num]
    
    # Add number of shared strains between these two species to the vector
    ani_list[length(ani_list)+1] <- pairs_matrix[sister1, sister2]
  }
  # Get mean ANI of all pairs, add it to a vector
  mean_ani_list[length(mean_ani_list)+1] <- mean(ani_list)
}

# Compute mean ANI for the real sister pairs
ani_list <- c()
for (genus in unique(df$Genus)){
  species = unique(df$Species_binomial[df$Genus == genus])
  ani <- pairs_matrix[species[1], species[2]]
  ani_list[length(ani_list)+1] <- ani
}

# add empirical mean number of shared strains to the vector with the fake ones
mean_ani_list[length(mean_ani_list)+1] <- 
  mean(ani_list)

# compute p-value (proportion of the reshuffled pairs that have mean strain sharing greater than or equal to the empirical one)e)
pvalue <- length(mean_ani_list[mean_ani_list >= mean(ani_list)])/length(mean_ani_list)

# plot the distribution
if(nuwt){
  nuwt_approach <- "(infection-conservative)"
} else {
  nuwt_approach <- "(NUWT-conservative)"
}

plotname <- paste0("histogram_mean_ANI", nuwt_approach)
for (plotnum in 1:3){
  if(plotnum==1){
    pdf(paste0(plotname, ".pdf"))
  } else if (plotnum==2){
    png(paste0(plotname, ".png"))
  } else {
    svg(paste0(plotname, ".svg"))
  }
  hist(mean_ani_list, xlim = c(min(mean_ani_list), max(mean_ani_list)),
       xlab = "Mean ANI between strains in the species pairs",
       main = paste0("Histogram of mean ANI between strains in species pairs\n",
                     nuwt_approach))
  abline(v = mean(ani_list), col = "red")
  if(nuwt){
    text(98, 3000, paste0("p-value = ", pvalue))
  } else {
    text(97.7, 3500, paste0("p-value = ", pvalue))
  }
  dev.off()
}
