# Supplementary Materials and code for Toro-Delgado _et al._ (in prep.)
This repository contains the Supplementary materials for the final project of the Master in Bioinformatics for Health Sciences of Pompeu Fabra University, by the student Eric Toro-Delgado. Refer to the main text for the relevant references to the supplementary figures and tables.

The contents of the `R_scripts` folder are as follows:
- `ANI_dxy_correlation.R`: script to conduct the correlation between the mean average nucleotide identity (ANI) of _Wolbachia_ strains found in a sister pair and the d_{xy} of between the host species of the pair.
- `make_species_pairs_ANI_matrix.R`: script to conduct the randomization test of whether _Wolbachia_ strains found in the species pairs have higher or lower ANI than expected by chance from randomly assigned pairs of hosts.
- `num_shared_strains_randomization.R`: script to conduct the randomization test of whether pairs of sister taxa share more _Wolbachia_ strains then expected by chance from randomly assigned pairs of hosts.
- `mcmcglmm_wolbachia_patterns_part1_code_only.R`: script to fit the Bayesian generalized linear model explaining the probability of sharing strains as a function of split time (in generations) and degree of symnpatry.
- `mcmcglmm_wolbachia_patterns_part1_split_time_mya_code_only.R`: the same as `mcmcglmm_wolbachia_patterns_part1_code_only.R`, but with split time in million years ago.
- `mcmcglmm_strain_number_code_only.R`: script to fit the Bayesian generalized linear model explaining the number of strains found in a species as a function of voltinism, wing index and number of flight months.
