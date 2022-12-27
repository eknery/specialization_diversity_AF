### packages
library(raster)
library(sp)
library(sf)
library(rgeos)
library(dismo)
library(ape)

#loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

# all species names
all_spp_names = sort(unique(spp_points$species))

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

# define states threshold
high_ths = 0.90
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "AF"
geo_states[af_percentage <= low_ths] = "other"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "AFother"
names(geo_states) = spp_count_domain$species

# sourcing other functions
source("2_sister_hypervolume/function_sister_pairs.R")

############################ sister area comparisons #######################

for (n in 1:length(phylo_trees)){
  one_tree = phylo_trees[[n]]
  # phylogenetic distance
  phylo_distance = cophenetic(one_tree)
  # sister taxa
  sister_taxa_list = sister_pairs(phylo_distance)
  ### sister divergence time
  sister_divergence = c()
  for(focal_sp in all_spp_names){
    focal_sp_dists= round(phylo_distance[focal_sp,],5)
    min_phylo_dist = round( min(phylo_distance[focal_sp,][phylo_distance[focal_sp,] != 0]), 5)
    sister_divergence = c(sister_divergence, min_phylo_dist)
  }
  ### sister areas
  sister_area_comparison = matrix(0, nrow=length(all_spp_names), ncol=4)
  for (i in 1:length(all_spp_names)){
    # focal species
    sp_name = all_spp_names[i]
    sp_points = spp_points[spp_points$species == sp_name, -1]
    sp_multi = st_multipoint(as.matrix(sp_points))
    sp_convex = sf::st_convex_hull(sp_multi)
    sister_area_comparison[i,1] = st_area(sp_convex)
    # sister taxa
    sister_name = sister_taxa_list[[sp_name]]
    sister_points = spp_points[spp_points$species %in% sister_name,-1]
    sister_multi = st_multipoint(as.matrix(sister_points))
    sister_convex = sf::st_convex_hull(sister_multi)
    sister_area_comparison[i,2] = st_area(sister_convex)
    # intersection area
    intersection = st_intersection(st_union(sp_convex),st_union(sister_convex))
    sister_area_comparison[i,3]  = st_area(intersection)
    # union area
    union = st_union(sp_convex, sister_convex)
    sister_area_comparison[i,4] = st_area(union)
  }
  # minimal area per sister pair
  min_area = apply(sister_area_comparison[,1:2], MARGIN = 1, FUN= min)
  # dataframe
  sister_area_comparison = data.frame(all_spp_names, sister_area_comparison, min_area, sister_divergence)
  colnames(sister_area_comparison) =c("species", "sp_area", "sister_area", "intersection", "union", "minimal_area","divergence_time")
  #exporting
  write.table(sister_area_comparison, paste("3_sister_geography/sister_area_comparisons/sister_area_comparison_", as.character(n), ".csv", sep="" ), sep=',', quote=F, row.names=F)
  # update me!
  print(paste("Time:", Sys.time(), "Loop iterarion:", as.character(n) ) )
}


######################## calculating sister RO metrics #########################

sister_ro_metrics = c()
for (n in 1:length(phylo_trees) ){ 
  sister_area_comparison = read.table(paste("3_sister_geography/sister_area_comparisons/sister_area_comparison_", as.character(n), ".csv", sep=""), sep=',', h=T)
  geo_groups = split(sister_area_comparison, geo_states)
  for (i in 1:length(geo_groups) ){
    group_name = names(geo_groups)[i]
    one_group = geo_groups[[i]]
    ro = one_group$intersection/one_group$minimal_area
    mean_ro = mean(ro)
    linear_model = lm(ro ~ one_group$divergence_time)
    intercept_ro =linear_model$coefficients["(Intercept)"]
    one_group_metric = c(group_name, mean_ro, intercept_ro)
    sister_ro_metrics = rbind(sister_ro_metrics, one_group_metric)
  }
}

sister_ro_metrics = data.frame(sister_ro_metrics)
colnames(sister_ro_metrics) = c("state", "mean_ro", "intercept_ro")
rownames(sister_ro_metrics) = NULL
sister_ro_metrics$mean_ro = as.numeric(sister_ro_metrics$mean_ro)
sister_ro_metrics$intercept_ro = as.numeric(sister_ro_metrics$intercept_ro)

#exporting
write.table(sister_ro_metrics, "3_sister_geography/sister_ro_metrics.csv", sep=',', quote=F, row.names=F)

############################## analyzing RO metrics ############################

sister_ro_metrics = read.table("3_sister_geography/sister_ro_metrics.csv", sep=',', h=T)

### summary by state
means = aggregate(sister_ro_metrics$intercept_ro, by = list(sister_ro_metrics$state), mean )
aggregate(sister_ro_metrics$intercept_ro, by = list(sister_ro_metrics$state), sd )

### testing differences
# source function
source("function_permutation_test.R")
# set factor
ro_metrics_split = split(sister_ro_metrics, f=sister_ro_metrics$state)
# loop over traits
all_pvalues = c()
for(i in 1:(length(ro_metrics_split)-1) ){
  for(j in (i+1):length(ro_metrics_split) ){
    # select data by state
    i_df = ro_metrics_split[[i]]
    j_df = ro_metrics_split[[j]]
    # state name
    i_name =unique(i_df$state)
    j_name =unique(j_df$state)
    # joing data into single df
    ij_df = rbind(i_df, j_df)
    pair_name = paste(i_name, j_name, sep="_")
    # set variable 
    response = ij_df$intercept_ro
    factor = ij_df$state
    # test
    pvalue = permutation_test(factor = factor, response = response, stats="mean", iter=999, out_dir="3_sister_geography",  name= paste(pair_name, ".tiff", sep=""))
    names(pvalue) = pair_name
    all_pvalues = c(all_pvalues, pvalue)
  }
}

all_pvalues
######################### geographic overlap ~ specialization #######################

spp_hvolumes = read.table("1_hypervolume_inference/spp_hvolumes.csv", sep=',', h=T)
hvolumes = spp_hvolumes$hvolume

sister_ro_metrics = c()
for (n in 1:length(phylo_trees) ){ 
  n = 1
  sister_area_comparison = read.table(paste("3_sister_geography/sister_area_comparisons/sister_area_comparison_", as.character(n), ".csv", sep=""), sep=',', h=T)
  ro = sister_area_comparison$intersection / sister_area_comparison$minimal_area
  linear_model = lm(ro ~ poly(hvolumes,2))
  
  plot(ro  ~ hvolumes )
  abline(linear_model)
}

hist( sqrt(ro) )
