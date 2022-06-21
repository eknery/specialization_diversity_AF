setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

### packages
library(raster)
library(rgeos)
library(hypervolume)
library(dismo)
library(ape)

#loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# loading geographic distribution
spp_geographic_distribution = read.table("0_data/spp_geographic_distribution.csv", sep=',', h=T)

### loading scale niche positions
hv_scale_positions = read.table("1_hypervolume_inference/hv_scale_positions.csv", sep=",", h=T)

# loading phylogenetic tree
mcc=read.tree("0_data/mcc_phylo.nwk")

# loading raster layers
ras1 = raster("0_data/rasters/temperature_diurnal_range")
ras2 = raster("0_data/rasters/precipitation_seasonality")
ras3 = raster("0_data/rasters/solar_radiation")
ras4 = raster("0_data/rasters/soil_pH.gri")
env_ras= stack(ras1,ras2,ras3,ras4)

# all species names
all_spp_names = sort(unique(spp_points$species))

# sourcing other functions
source("2_sister_hypervolume/function_sister_pairs.R")

################################## data preparation ###########################

### croping rasters
#set extent
ext=extent(-100, -33, -34, 24)
allras_crop=crop(env_ras[[1]], ext)
for (i in 2:length(env_ras@layers) ){
  ras_crop =  crop(env_ras[[i]], ext)
  allras_crop= raster::stack(allras_crop, ras_crop)
}
crop_env_ras = allras_crop


### keeping mean and sd values
ras_means= cellStats(crop_env_ras, stat='mean', na.rm=TRUE)
ras_sds = cellStats(crop_env_ras, stat='sd', na.rm=TRUE)

### converting rasters to z-scores
scale_env_ras = crop_env_ras[[c(1:nlayers(env_ras))]]
for (i in 1:nlayers(scale_env_ras)) { 
  scale_env_ras[[i]] = (scale_env_ras[[i]] - cellStats(scale_env_ras[[i]], 'mean')) / cellStats(scale_env_ras[[i]], 'sd') 
}

### treating spp occurrences
# extracting spp z-values
scale_spp_env = extract(scale_env_ras,spp_points[,2:3])
scale_spp_env = data.frame(spp_points$species, scale_spp_env)
colnames(scale_spp_env)[1] = "species"

# removing spp NAs
spp_nas = c()
for (i in 2:length(scale_spp_env[1,])){
  spp_nas= c(spp_nas, which(is.na(scale_spp_env[,i])) )}
if (length(spp_nas) > 0){
  scale_spp_env = scale_spp_env[-spp_nas,]}

############################ loading best model threshold ######################

# loading threshold values
mean_hv_performance=read.table("1_hypervolume_inference/mean_hv_performance.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# getting best thresholds for each species
best_ths = rep(NA, length(all_spp_names))
for(sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  mean_sp_performance = mean_hv_performance[mean_hv_performance$species==sp_name,]
  sp_th = mean_sp_performance$threshold[mean_sp_performance$tss == max(mean_sp_performance$tss)] 
  best_ths[index] = sp_th[1]
}
names(best_ths) = all_spp_names

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

################################# sister hv comparison ##########################

sister_hv_comparison = matrix(0, nrow=length(all_spp_names), ncol=4)
for (sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  sp_data = scale_spp_env[scale_spp_env$species== sp_name, 2:5]
  sp_ths = best_ths[names(best_ths) == sp_name]
  sp_band = estimate_bandwidth(sp_data)
  sp_hv = hypervolume_gaussian(sp_data,  samples.per.point = ceiling((10^(3 + sqrt(ncol(sp_data))))/nrow(sp_data)), 
                               kde.bandwidth = sp_band, sd.count = 3, chunk.size = 100,
                               quantile.requested = sp_ths, quantile.requested.type = "probability")
  sister_hv_comparison[index,1] = sp_hv@Volume
  sister_data = scale_spp_env[scale_spp_env$species %in% sister_taxa_list[[sp_name]], 2:5]
  sister_ths = median(best_ths[names(best_ths) %in% sister_taxa_list[[sp_name]]])
  sister_band = estimate_bandwidth(sister_data)
  sister_hv = hypervolume_gaussian(sister_data,  samples.per.point = ceiling((10^(3 + sqrt(ncol(sister_data))))/nrow(sister_data)), 
                               kde.bandwidth = sister_band, sd.count = 3, chunk.size = 100,
                               quantile.requested = sister_ths, quantile.requested.type = "probability")
  sister_hv_comparison[index,2] = sister_hv@Volume
  # getting intersection and union
  hv_set = hypervolume_set(sp_hv, sister_hv, check.memory = F)
  sister_hv_comparison[index,3] = hv_set[[3]]@Volume # taking intersection
  sister_hv_comparison[index,4] = hv_set[[4]]@Volume # taking union
}

# minimal hv per sister pair
min_hv = apply(sister_hv_comparison[,1:2], MARGIN = 1, FUN= min)

# dataframe
sister_hv_comparison = data.frame(spp_geographic_distribution$state, all_spp_names, sister_hv_comparison, min_hv, sister_divergence_time$sister_divergence)
colnames(sister_hv_comparison) =c("state", "species", "sp_hvolume", "sister_hvolume", "intersection", "union", "minimal_hvolume", "divergence_time")

# exporting
write.table(sister_hv_comparison, "2_sister_hypervolume/sister_hv_comparison.csv", sep=",", quote=F, row.names = F)

#################################### sister hv position #########################

distance_to_sister = c()
for (i in 1:length(all_spp_names)){
  # focal sp position
  sp_name = all_spp_names[i]
  sp_positions = hv_scale_positions[hv_scale_positions$species == sp_name,-1]
  # sister sp position
  sister_name = sister_taxa_list[[sp_name]]
  sister_positions = hv_scale_positions[hv_scale_positions$species %in% sister_name,-1]
    if( nrow(sister_positions) != 1 ){ # take median if more than one sister
      sister_positions = apply(sister_positions, MARGIN=2, FUN=median)
    }
  # eucledian distance
  distance_to_sister[i] = sqrt(sum((sp_positions - sister_positions)^2))
}

sister_hv_distance = data.frame(spp_geographic_distribution$state, all_spp_names, distance_to_sister, sister_divergence_time$sister_divergence)
colnames(sister_hv_distance) = c("state", "species", "distance_to_sister", "divergence_time")

#exporting
write.table(sister_hv_distance, "2_sister_hypervolume/sister_hv_distance.csv", sep=',', quote=F, row.names=T)
