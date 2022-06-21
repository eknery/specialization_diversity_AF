setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

#loading spp coordinates
phylo_distance=read.table("phylo_distance.csv", header =T, sep=",",  na.strings = "NA", fill=T)
str(phylo_distance)

#loading phylogenetic tree
library(ape)
mcc=read.tree("mcc_phylo.nwk")

### list of sister taxa
all_spp_names = rownames(phylo_distance)

sister_taxa_list = vector('list', length(all_spp_names))
names(sister_taxa_list) = all_spp_names

for(focal_sp in all_spp_names){
  focal_sp_dists= round(phylo_distance[focal_sp,],5)
  min_phylo_dist = round( min(phylo_distance[focal_sp,][phylo_distance[focal_sp,] != 0]), 5)
  sister_index = which(focal_sp_dists == min_phylo_dist)
  sister_taxa_list[[focal_sp]] = colnames(focal_sp_dists)[sister_index]
}




