############################# calculating phylogenetic distance #######################

#loading package
library(ape)

# loading phylogenetic tree
mcc=read.tree("0_data/mcc_phylo.nwk")

# species names
all_spp_names = unique(mcc$tip.label)
all_spp_names =all_spp_names[order(all_spp_names)]

init_phylo_distance = cophenetic(mcc)

row_phylo_distance = init_phylo_distance[which(rownames(init_phylo_distance) == all_spp_names[1]),]
for (i in 2:length(all_spp_names)){
  next_row = init_phylo_distance[which(rownames(init_phylo_distance) == all_spp_names[i]),]
  row_phylo_distance= rbind(row_phylo_distance, next_row)
}
rownames(row_phylo_distance) = all_spp_names

col_phylo_distance = row_phylo_distance[ ,which(colnames(row_phylo_distance) == all_spp_names[1])]
for (i in 2:length(all_spp_names)){
  next_col = row_phylo_distance[ ,which(colnames(row_phylo_distance) == all_spp_names[i])]
  col_phylo_distance = cbind(col_phylo_distance, next_col)
}
colnames(col_phylo_distance) = all_spp_names

# organizing
phylo_distance = col_phylo_distance

# exporting
write.table(phylo_distance, "0_data/phylo_distance.csv", sep=",", quote=F, row.names = T)

