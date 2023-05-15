#!/usr/bin/env Rscript

# Script to test if there is a correlation between d_xy and ANI

# The rationale is that if wolbachia tends to covdiverge with its host, there should
# be a correlation between divergence of sister host taxa and divergence of the strains
# that are codiverging with these hosts (that is, if they have not been 
# replaced by another strain)


#### LIBS ####
library(tidyverse)

#### RAW DATA ####
# since we have two version of the data, we run it for the infection-conservative version
# and then for the NUWT-conservative version

nuwt <- F

# Samples list
if(nuwt){
  df <- read.csv("../Samples_all_wolbachia_presence_absence_infection_conservative.csv")
} else {
  df <- read.csv("../Samples_all_wolbachia_presence_absence_NUWT_conservative.csv")
}

# Data with the species pairs traits from Ebdon et al. (2020)
if(nuwt){
  pairs_traits <- read.csv("../Species_pairs_traits_infection_conservative.csv")
} else {
  pairs_traits <- read.csv("../Species_pairs_traits_NUWT_conservative.csv")
}

# ANI values for all genomes computed with fastANI
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

# Keep only pairs where both are infected (otherwise we can't get an ANI)
pairs_traits <- pairs_traits[pairs_traits$Both.infected == 1,]

# Keep only species pairs where we have the traits
pairs_traits <- pairs_traits[!is.na(pairs_traits$d_xy),]

# Clean up genus names
pairs_traits$Genus <- gsub("*", "", pairs_traits$Genus, fixed = T)

# Expand genus names in Species2
for (rownum in 1:nrow(pairs_traits)){
  pairs_traits$Species.2[rownum] <- gsub("[A-Z]\\.", pairs_traits$Genus[rownum], pairs_traits$Species.2[rownum])
}


# Build mean ANI matrix of N x N
species_list = sort(c(pairs_traits$Species.1, pairs_traits$Species.2))
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

# add mean ANI of each pair to the pairs df
for (rownum in 1:nrow(pairs_traits)){
  pairs_traits$ANI[rownum] <- pairs_matrix[pairs_traits$Species.1[rownum], pairs_traits$Species.2[rownum]]
}

# build correlation of d_xy and ANI
if(nuwt){
  nuwt_approach <- "(infection-conservative)"
} else {
  nuwt_approach <- "(NUWT-conservative)"
}

for (plotnum in 1:3){
  if (plotnum == 1){
    pdf(paste0("ANI_dxy_correlation", nuwt_approach, ".pdf"))
  } else if (plotnum==2){
    png(paste0("ANI_dxy_correlation", nuwt_approach, ".png"))
  } else {
    svg(paste0("ANI_dxy_correlation", nuwt_approach, ".svg"))
  }
  
  plot(pairs_traits$ANI~pairs_traits$d_xy,
       ylab = paste0("Wolbachia", " average nucleotide identity (ANI)"),
       xlab = "Host d_xy",
       main = paste0("Correlation between host and\n", "Wolbachia", " divergence ", nuwt_approach)
  )
  
  model <- lm(pairs_traits$ANI~pairs_traits$d_xy)
  modelsum <- summary(model)
  
  corrtest <- cor.test(pairs_traits$ANI, pairs_traits$d_xy, method = "pearson")
  abline(a = modelsum$coefficients[1, 1], b = modelsum$coefficients[2, 1])
  text(0.07, 99.5, paste0("R2 = ", as.character(round(corrtest$estimate, digits = 2)), "; p-value = ", as.character(round(corrtest$p.value, digits = 3))))
  dev.off()
}
