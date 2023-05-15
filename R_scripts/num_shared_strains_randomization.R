#!/usr/bin/env Rscript

# Script to test if the number of strains shared by species pairs is higher or 
# lower than expected by chance.

# The rationale is that if Wolbachia is more easily spread across more closely 
# related species, then these species pairs will share strains more often that
# random species pairs, but if differing Wolbachia infection status favours speciation
# at early stages, then the young sister pairs will tend to share Wolbachia strains
# less often than random pairs

# all this makes sense for pairs where at least one has wolbachia, otherwise
# we can't ask if they share strains or not if none are infected

#### LIBRARIES ####
library(tidyverse)

# RAW DATA
# since we have samples of uncertain infection status, we run this for two versions
# of the data, one "infection-conservative" and one "NUWT-conservative"

nuwt = F

if (nuwt){
  df <- read.csv("../Samples_all_wolbachia_presence_absence_infection_conservative.csv")
} else {
  df <- read.csv("../Samples_all_wolbachia_presence_absence_NUWT_conservative.csv")
}




#### MAIN ####
# Add species binomial name
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

df_all <- df
rm(df)

# Make dataframe without considering pairs where one or both have no Wolbachia
df_both_wolbachia <- df_all[df_all$Genus %in% c("Colias", "Erebia", "Gonepteryx", 
                                        "Iphiclides", "Lasiomamta", "Melanargia", 
                                        "Thymelicus"),]

# Make dataframe considering pairs where at least one of them has some Wolbachia
df_one_or_more_wolbachia <- df_all[df_all$Genus %in% c("Colias", "Erebia", "Gonepteryx", 
                                               "Iphiclides", "Lasiomamta", 
                                               "Melanargia", "Thymelicus", "Pyrgus",
                                               "Pontia", "Polyommatus"),]

# run analysis for all three dataframes
dataframes_list <- list(df_all, df_one_or_more_wolbachia, df_both_wolbachia)

if (nuwt){
  nuwt_approach = "(infection-conservative, "
} else {
  nuwt_approach = "(NUWT-conservative, "
}

for (dfnum in 1:length(dataframes_list)){
  df <- dataframes_list[[dfnum]]
  
  pairs_considered_approach <- ifelse(dfnum==1, "all pairs)", 
                                      ifelse(dfnum==2, "at least one species infected)",
                                             "both species infected)"))
  
  xcoord <- ifelse(dfnum==1, 0.3, ifelse(dfnum==2, 0.6, 0.9))
  ycoord <- ifelse(dfnum==1, 12500, ifelse(dfnum==2, 10000, 8000))
  
  # Build matrix of N x N
  species_list = unique(df$Species_binomial)
  species_num = length(species_list)
  pairs_matrix <- as.data.frame(matrix(nrow = species_num, ncol = species_num))
  colnames(pairs_matrix) = species_list
  rownames(pairs_matrix) = species_list
  
  for (species1 in species_list){
    for (species2 in species_list){
      strains_sp1 = unique(df$Closest.reference.Wolbachia[df$Species_binomial == species1 & df$Wolbachia.presence == "Y"])
      strains_sp2 = unique(df$Closest.reference.Wolbachia[df$Species_binomial == species2 & df$Wolbachia.presence == "Y"])
      common_strains <- intersect(strains_sp1, strains_sp2)
      common_strains <- na.omit(common_strains)
      common_strains_num <- length(common_strains)
      pairs_matrix[species1, species2] <- common_strains_num
    }
  }
  
  # Build fake sets of sister pairs and calculate average number of strains shared in each pair
  # Vector to store mean number of strains shared in each fake dataset
  mean_num_shared_strains_list <- c()
  
  for (replicate in 1:19999){
    # Create vector to store number of shared strains in each of the fake pairs
    num_shared_strains_list <- c()
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
      num_shared_strains_list[length(num_shared_strains_list)+1] <- pairs_matrix[sister1, sister2]
    }
    # Get mean number of strains shared in the pairs, add it to a vector
    mean_num_shared_strains_list[length(mean_num_shared_strains_list)+1] <- mean(num_shared_strains_list)
  }
  
  # Compute mean number of strains shared for the real sister pairs
  num_shared_strains_list <- c()
  for (genus in unique(df$Genus)){
    species = unique(df$Species_binomial[df$Genus == genus])
    num_shared_strains <- pairs_matrix[species[1], species[2]]
    num_shared_strains_list[length(num_shared_strains_list)+1] <- num_shared_strains
  }
  
  # add empirical mean number of shared strains to the vector with the fake ones
  mean_num_shared_strains_list[length(mean_num_shared_strains_list)+1] <- 
    mean(num_shared_strains_list)
  
  # compute p-value (proportion of the reshuffled pairs that have mean strain sharing greater than or equal to the empirical one)e)
  pvalue <- length(mean_num_shared_strains_list[mean_num_shared_strains_list >= mean(num_shared_strains_list)])/length(mean_num_shared_strains_list)
  
  # plot the distribution
  plot_name <- gsub(",", "", 
                    gsub(" ", "_", 
                         paste0("Histogram_mean_number_shared_strains_", nuwt_approach, pairs_considered_approach)
                    )
  )
  for (plotnum in 1:3) {
    if(plotnum == 1){
      pdf(paste0(plot_name, ".pdf"))
    } else if (plotnum == 2) {
      svg(paste0(plot_name, ".svg"))
    } else {
      png(paste0(plot_name, ".png"))
    }
    hist(mean_num_shared_strains_list, 
         xlim = c(0, max(mean_num_shared_strains_list)),
         xlab = "Mean number of shared strains per species pair",
         main = paste0("Histogram of mean number of shared strains per species pair\n", nuwt_approach, pairs_considered_approach))
    abline(v = mean(num_shared_strains_list), col = "red")
    text(xcoord, ycoord, labels = paste0("p-value = ", pvalue))
    dev.off()
  }
}
