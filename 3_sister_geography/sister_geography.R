setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

### packages
library(raster)
library(sp)
library(sf)
library(rgeos)
library(dismo)
library(ape)

#loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# loading geographic distribution
spp_geographic_distribution = read.table("0_data/spp_geographic_distribution.csv", sep=',', h=T)

# loading phylogenetic tree
mcc=read.tree("0_data/mcc_phylo.nwk")

# all species names
all_spp_names = sort(unique(spp_points$species))

# altitude raster
altitude_ras = raster("0_data/rasters/altitude")

# sourcing other functions
source("2_sister_hypervolume/function_sister_pairs.R")

############################# sister taxa and divergence time ############################

# phylogenetic distance
phylo_distance = cophenetic(mcc)

# sister taxa
sister_taxa_list = sister_pairs(phylo_distance)

### sister divergence time
sister_divergence = c()
for(focal_sp in all_spp_names){
  focal_sp_dists= round(phylo_distance[focal_sp,],5)
  min_phylo_dist = round( min(phylo_distance[focal_sp,][phylo_distance[focal_sp,] != 0]), 5)
  sister_divergence = c(sister_divergence, min_phylo_dist)
}

sister_divergence_time = data.frame(all_spp_names,sister_divergence)
colnames(sister_divergence_time)[1] = "species"

############################ sister area comparison #######################

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
sister_area_comparison = data.frame(spp_geographic_distribution$state, all_spp_names, sister_area_comparison, min_area, sister_divergence_time$sister_divergence)
colnames(sister_area_comparison) =c("state", "species", "sp_area", "sister_area", "intersection", "union", "minimal_area","divergence_time")

#exporting
write.table(sister_area_comparison, "3_sister_geography/sister_area_comparison.csv", sep=',', quote=F, row.names=T)

######################## sister geo distance ###############################

sister_geo_distance = c()
for (i in 1:length(all_spp_names)){
  # focal species
  sp_name = all_spp_names[i]
  sp_points = spp_points[spp_points$species == sp_name, -1]
  sp_sf = st_as_sf(sp_points, coords=c(1,2))
  # sister taxa
  sister_name = sister_taxa_list[[sp_name]]
  sister_points = spp_points[spp_points$species %in% sister_name,-1]
  sister_sf = st_as_sf(sister_points, coords=c(1,2))
  # calculating distance
  sister_geo_distance[i] = median(st_distance(sp_sf, sister_sf))
}

# dataframe
sister_geo_distance = data.frame(spp_geographic_distribution$state, all_spp_names, sister_geo_distance, sister_divergence_time$sister_divergence)
colnames(sister_geo_distance) = c("state", "species", "distance_to_sister", "divergence_time")

#exporting
write.table(sister_geo_distance, "3_sister_geography/sister_geo_distance.csv", sep=',', quote=F, row.names=T)

