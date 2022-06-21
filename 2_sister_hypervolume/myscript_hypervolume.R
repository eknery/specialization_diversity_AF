setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

### packages
library(raster)
library(sp)
library(sf)
library(spatialEco)
library(rgeos)
library(hypervolume)
library(dismo)
library(ape)


#loading spp coordinates
spp_points=read.table("0_data/spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)

# loading phylogenetic tree
mcc=read.tree("0_data/mcc_phylo.nwk")

# loading raster layers
ras1 = raster("0_data/rasters/temperature_diurnal_range")
ras2 = raster("0_data/rasters/precipitation_seasonality")
ras3 = raster("0_data/rasters/solar_radiation")
ras4 = raster("0_data/rasters/soil_pH.gri")
env_ras= stack(ras1,ras2,ras3,ras4)

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


############################# calculating phylogenetic distance #######################

init_phylo_distance = cophenetic(mcc)

row_phylo_distance = init_phylo_distance[which(rownames(init_phylo_distance) == rownames(hv_similarity)[1]),]
for (i in 2:nrow(hv_similarity)){
  next_row = init_phylo_distance[which(rownames(init_phylo_distance) == rownames(hv_similarity)[i]),]
  row_phylo_distance= rbind(row_phylo_distance, next_row)
}
rownames(row_phylo_distance) = rownames(hv_similarity)

col_phylo_distance = row_phylo_distance[ ,which(colnames(row_phylo_distance) == colnames(hv_similarity)[1])]
for (i in 2:ncol(hv_similarity)){
  next_col = row_phylo_distance[ ,which(colnames(row_phylo_distance) == colnames(hv_similarity)[i])]
  col_phylo_distance = cbind(col_phylo_distance, next_col)
}
colnames(col_phylo_distance) = colnames(hv_similarity)
phylo_distance = col_phylo_distance

write.table(phylo_distance, "phylo_distance.csv", sep=",", quote=F, row.names = T)

############################## sister species comparisons ###########################

### sister taxa and divergence time
sister_taxa_list = vector('list', length(all_spp_names))
sister_divergence = c()
names(sister_taxa_list) = all_spp_names
for(focal_sp in all_spp_names){
  focal_sp_dists= round(phylo_distance[focal_sp,],5)
  min_phylo_dist = round( min(phylo_distance[focal_sp,][phylo_distance[focal_sp,] != 0]), 5)
  sister_divergence = c(sister_divergence, min_phylo_dist)
  sister_index = which(focal_sp_dists == min_phylo_dist)
  sister_taxa_list[[focal_sp]] = colnames(focal_sp_dists)[sister_index]
}

sister_divergence_time = data.frame(all_spp_names,sister_divergence)
colnames(sister_divergence_time)[1] = "species"

# exporting
write.table(sister_divergence_time, "sister_divergence_time.csv", sep=",", quote=F, row.names = F)

### hv similarity
sister_hv_similarity = rep(NA, length(all_spp_names))
for (sp_name in all_spp_names){
  index = which(sp_name == all_spp_names)
  sp_data = scale_spp_env[scale_spp_env$species== sp_name, 2:5]
  sp_ths = best_ths[names(best_ths) == sp_name]
  sp_band = estimate_bandwidth(sp_data)
  sp_hv = hypervolume_gaussian(sp_data,  samples.per.point = ceiling((10^(3 + sqrt(ncol(sp_data))))/nrow(sp_data)), 
                               kde.bandwidth = sp_band, sd.count = 3, chunk.size = 100,
                               quantile.requested = sp_ths, quantile.requested.type = "probability")
  sister_data = scale_spp_env[scale_spp_env$species %in% sister_taxa_list[[sp_name]], 2:5]
  sister_ths = median(best_ths[names(best_ths) %in% sister_taxa_list[[sp_name]]])
  sister_band = estimate_bandwidth(sister_data)
  sister_hv = hypervolume_gaussian(sister_data,  samples.per.point = ceiling((10^(3 + sqrt(ncol(sister_data))))/nrow(sister_data)), 
                               kde.bandwidth = sister_band, sd.count = 3, chunk.size = 100,
                               quantile.requested = sister_ths, quantile.requested.type = "probability")
  hv_set = hypervolume_set(sp_hv, sister_hv, check.memory = F)
  sister_hv_similarity[index] = as.numeric(hypervolume_overlap_statistics(hv_set)[1])
}
sister_hv_similarity = data.frame(all_spp_names, sister_hv_similarity)
colnames(sister_hv_similarity) =c("species", "similirarity_to_sister")

### niche dissimilarity
sister_hv_dissimilarity = sister_hv_similarity
sister_hv_dissimilarity[,2] = 1- sister_hv_dissimilarity[,2]
colnames(sister_hv_dissimilarity)[2] = "dissimilarity_to_sister"

# exporting
write.table(sister_hv_similarity, "sister_hv_similarity.csv", sep=",", quote=F, row.names = F)

### geographic distance
sister_geo_distance = c()
for (sp_name in all_spp_names){
  sp_focal = spp_points[spp_points$species == sp_name, 2:3]
  sister = spp_points[spp_points$species %in% sister_taxa_list[[sp_name]],2:3]
  med_point_distance = median(pointDistance(sp_focal, sister, lonlat=T)*10^-3)
  sister_geo_distance = rbind(sister_geo_distance, c(sp_name, med_point_distance))
}
sister_geo_distance = data.frame(sister_geo_distance)
colnames(sister_geo_distance) =c("species", "distance_to_sister")

#exporting
write.table(sister_geo_distance, "sister_geo_distance.csv", sep=",", quote=F, row.names = F)

### distribution vector
distribution = sister_hv_similarity$species
distribution[sister_hv_similarity$species %in% AF] = "AF-endemic"
distribution[!sister_hv_similarity$species %in% AF] = "widespread"

### comparing

library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)

sister_data = data.frame(distribution, sister_hv_dissimilarity, as.numeric(sister_geo_distance[,2]), sister_divergence_time[,2] )
colnames(sister_data)[4:5] = c("geodistance_to_sister", "divergence")

ggplot(data= sister_data, aes(x=distribution, y=dissimilarity_to_sister/divergence, fill=distribution)) +
  geom_point(aes(y=dissimilarity_to_sister/divergence, color=distribution), position = position_jitter(width = 0.07), size = 2, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")

ggplot(data= sister_data, aes(x=distribution, y=geodistance_to_sister/divergence, fill=distribution)) +
  geom_point(aes(y=geodistance_to_sister/divergence, color=distribution), position = position_jitter(width = 0.07), size = 2, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  ylim(c(0,1500))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")


############################ overall species comparisons ############################

### setting function to convert matrix to dataframe
matrix_to_df = function(matrix){
  # AF matrix
  af_matrix = matrix[rownames(matrix) %in% AF, colnames(matrix) %in% AF]
  af_values = c()
  for (i in 1:nrow(af_matrix)){
    for(j in i:ncol(af_matrix)){
      if(i == j){next}
      af_values= c(af_values, af_matrix[i,j])
    }
  }
  
  # WS matrix
  ws_matrix = matrix[!rownames(matrix) %in% AF, !colnames(matrix) %in% AF]
  ws_values = c()
  for (i in 1:nrow(ws_matrix)){
    for(j in i:ncol(ws_matrix)){
      if(i == j){next}
      ws_values= c(ws_values, ws_matrix[i,j])
    }
  }
  
  # organizing into dataframe
  af_name = rep("AF-endemic", length(af_values))
  ws_name = rep("widespread", length(ws_values))
  distribution = c(af_name,ws_name)
  all_values = c(af_values,ws_values)
  df = data.frame(distribution, all_values)
  return(df)
}

### niche similarity
hv_similarity = read.table("hv_similarity.csv", sep=",", h=T)
hv_similarity_df = matrix_to_df(hv_similarity)

# hv dissimilarity
hv_dissimilarity_df = hv_similarity_df
hv_dissimilarity_df[,2] = 1 - hv_dissimilarity_df[,2] 
colnames(hv_dissimilarity_df)[2] = "Jaccard_dissimilarity"

ggplot(data= hv_dissimilarity_df, aes(x=distribution, y=Jaccard_dissimilarity, fill=distribution)) +
  geom_point(aes(y=Jaccard_dissimilarity, color=distribution), position = position_jitter(width = 0.07), size = 1, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")


### median point distance
med_point_distance = read.table("med_point_distance.csv", sep=",", h=T)
med_point_distance_df = matrix_to_df(med_point_distance)
colnames(med_point_distance_df)[2] = "geographic_distance"

ggplot(data= med_point_distance_df,aes(x=distribution, y=geographic_distance, fill=distribution)) +
  geom_point(aes(y=geographic_distance, color=distribution), position = position_jitter(width = 0.07), size = 1, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  labs(y="pairwise geographic distance (km)", x="distribution")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")


### phylogenetic distance
phylo_distance = read.table("phylo_distance.csv", sep=",", h=T)
phylo_distance_df = matrix_to_df(phylo_distance)
colnames(phylo_distance_df)[2] = "time_of_divergence"

ggplot(data= phylo_distance_df,aes(x=distribution, y=time_of_divergence, fill=distribution)) +
  geom_point(aes(y=time_of_divergence, color=distribution), position = position_jitter(width = 0.07), size = 1, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")


### all pariwise comparisons
pair_df = data.frame(hv_dissimilarity_df, med_point_distance_df[,2], phylo_distance_df[,2])
colnames(pair_df)[3:4] = c("geographic_distance", "time_divergence")

# boxplot
ggplot(data= pair_df,aes(x=distribution, y=Jaccard_dissimilarity/time_divergence , fill=distribution)) +
  geom_point(aes(y=Jaccard_dissimilarity/time_divergence, color=distribution), position = position_jitter(width = 0.07), size = 1, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  ylim(c(0,0.5)) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")

### continuous plots
ggplot(data= pair_df, aes(y=Jaccard_dissimilarity/time_divergence, x=geographic_distance/time_divergence, fill=distribution)) +
  geom_point(aes(color=distribution), size = 2, alpha = 0.25) +
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  xlim(c(0,300))+
  ylim(c(0,0.3))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")

ggplot(data= pair_df, aes(x=time_divergence, y=geographic_distance, fill=distribution)) +
    geom_point(aes(color=distribution), size = 2, alpha = 0.25) +
    scale_fill_manual(values=c("green3","magenta"))+
    scale_colour_manual(values=c("green3","magenta"))+
    theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")
  
############################ EXTRA comparisons ###############################

#################### calculating niche position distance #######################

### loading scale niche positions
hv_scale_positions = read.table("hv_scale_positions.csv", sep=",", h=T)
positions= hv_scale_positions[,2:5]
rownames(positions) = hv_scale_positions$species

position_distance = as.matrix(dist(positions, "euclidean"))

#exporting
write.table(position_distance, "position_distance.csv", sep=',', quote=F, row.names=T)

############################ calculating area overlap #######################

### setting result matrix
area_similarity = data.frame( matrix(NA, nrow=length(all_spp_names), ncol=length(all_spp_names) ) )
colnames(area_similarity) = all_spp_names
rownames(area_similarity) = all_spp_names

for (i in 1:length(all_spp_names)){
  sp_name_i = all_spp_names[i]
  sp_points_i = spp_points[spp_points$species == sp_name_i, 2:3]
  st_i = st_multipoint(as.matrix(sp_points_i))
  convex_i = sf::st_convex_hull(st_i)
  for (j in i:length(all_spp_names)){
    sp_name_j = all_spp_names[j]
    sp_points_j = spp_points[spp_points$species == sp_name_j, 2:3]
    st_j = st_multipoint(as.matrix(sp_points_j))
    convex_j = sf::st_convex_hull(st_j)
    # intersection area
    inter_ij = st_intersection(st_union(convex_i),st_union(convex_j))
    inter_area_ij = st_area(inter_ij)
    # union area
    union_ij = st_union(convex_i, convex_j)
    union_area_ij = st_area(union_ij)
    # keeping proportional area
    area_similarity[i,j] = inter_area_ij/union_area_ij
  }
}

#exporting
write.table(area_similarity, "area_similarity.csv", sep=',', quote=F, row.names=T)

############################## ALTITUDE ####################################

spp_altitude = data.frame(spp_points$species, raster::extract(alt_ras, spp_points[,2:3]) )
med_spp_alt = aggregate(spp_altitude[,2], by=list(spp_altitude[,1]), function(x){median(x, na.rm=T)})

distribution = med_spp_alt$Group.1
distribution[distribution %in% AF] = "AF-endemic"
distribution[distribution!= "AF-endemic"] = "widespread"

med_spp_alt = data.frame(distribution, med_spp_alt)
colnames(med_spp_alt) = c("distribution", "species", "altitude")

altitude = med_spp_alt$altitude
names(altitude) = med_spp_alt$species

############################# niche structure ##############################


##### correlation structure

scale_env_values = scale_spp_env[,2:5]

for (i in ncol(scale_env_values)){
cor.test(scale_env_values[,1], scale_env_values[,4], method= "spearman")
}

###### PCA 
correlation = cor(scale_spp_env[,2:5])
pca = prcomp(scale_spp_env[,2:5], center = T)
stdev = pca$sdev / sum(pca$sdev)
pc_number =c("1", "2", "3", "4")
stdev = data.frame(pc_number, stdev)

### boots values
boot_stdev_df = data.frame(matrix(0, nrow=1, ncol=4))
boot_rotation_list = vector('list', 99)

for (i in 1:99){
  boot_num = sample(x=1:nrow(scale_spp_env[,2:5]), size = nrow(scale_spp_env[,2:5]), replace = T)
  boot_env_values = scale_spp_env[boot_num,2:5]
  boot_pca = prcomp(boot_env_values, center = T)
  boot_rotation_list[[i]] = boot_pca$rotation
  boot_stdev = boot_pca$sdev/ sum(boot_pca$sdev)
  boot_stdev_df = rbind(boot_stdev_df, boot_stdev)
}
boot_stdev_df = boot_stdev_df[-1,]

# generating stdev intervals
boot_stdev_means = apply(boot_stdev_df, 2,FUN=mean)
boot_stdev_quantiles = apply(boot_stdev_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
boot_stdev_quantiles= t(boot_stdev_quantiles)
stdev = data.frame(stdev, boot_stdev_quantiles)

###### randomization procedure 
rand_correlation_list = vector('list', 99)
rand_stdev_df = data.frame(matrix(0, nrow=1, ncol=4))
rand_rotation_list = vector('list', 99)

for (i in 1:99){
  rand_env_values = scale_spp_env[,2:5]
  for (j in ncol(rand_env_values) ){
    rand_num = sample(x=1:nrow(rand_env_values), replace=F)
    rand_env_values[,j] = rand_env_values[rand_num, j]
  }
  rand_correlation_list[[i]] = cor(rand_env_values)
  rand_pca = prcomp(rand_env_values, center = T)
  rand_rotation_list[[i]] = rand_pca$rotation
  rand_stdev = rand_pca$sdev/ sum(rand_pca$sdev)
  rand_stdev_df = rbind(rand_stdev_df, rand_stdev)
}
rand_stdev_df = rand_stdev_df[-1,]

### collecting stdev statistics
rand_stdev_means = apply(rand_stdev_df, 2,FUN=mean)
rand_stdev_quantiles = apply(rand_stdev_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
rand_stdev_stats= t(rbind(rand_stdev_means, rand_stdev_quantiles))
rand_stdev_stats = data.frame(pc_number, rand_stdev_stats)

### plotting
tiff("PC_variance.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data=rand_stdev_stats, aes(x=pc_number, y=rand_stdev_means*100) ) +
  geom_point(data=rand_stdev_stats, aes(y=rand_stdev_means*100), color="gray", size = 2) +
  geom_errorbar(aes(ymin=X5.*100, ymax=X95.*100), col="gray", size=1, width=0.25)+
  geom_point(data=stdev, aes(y=stdev*100), color="red", size = 1.5) +
  geom_errorbar(data=stdev, aes(ymin=X5.*100, ymax=X95.*100), col="red", size=1, width=0)+
  labs(x="PC axis", y="% explained variance")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12))
dev.off()

### collecting rotation statistics
rand_rotation_df = rand_rotation_list[[1]][,2]
for (i in 2:length(rand_rotation_list) ){
  rand_rotation_df = rbind(rand_rotation_df , rand_rotation_list[[i]][,2])
}
rand_rotation_means = apply(rand_rotation_df, 2,FUN=mean)
rand_rotation_quantiles = apply(rand_rotation_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
rand_rotation_stats= t(rbind(rand_rotation_means, rand_rotation_quantiles))
variable = c("1","2","3","4")
rand_rotation_stats = data.frame(variable, rand_rotation_stats)

# observed rotation
boot_rotation_df = boot_rotation_list[[1]][,2]
for (i in 2:length(boot_rotation_list) ){
  boot_rotation_df = rbind(boot_rotation_df , boot_rotation_list[[i]][,2])
}

boot_rotation_quantiles = apply(boot_rotation_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
boot_rotation_stats= t(rbind(pca$rotation[,2], boot_rotation_quantiles))
variable = c("1","2","3","4")
boot_rotation_stats = data.frame(variable, boot_rotation_stats)
colnames(boot_rotation_stats)[2] = "rotation"


### ploting
tiff("PC2_loadings.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data=rand_rotation_stats, aes(x=variable, y=rand_rotation_means) ) +
  geom_point(data=rand_rotation_stats, aes(y=rand_rotation_means), color="gray", size = 2) +
  geom_errorbar(aes(ymin= X5., ymax= X95.), col="gray", size=1, width=0.25)+
  geom_point(data=boot_rotation_stats, aes(y=rotation), color="red", size = 1.5) +
  geom_errorbar(data=boot_rotation_stats, aes(ymin= X5., ymax= X95.), col="red", size=1, width=0)+
  labs(x="variable", y="PC2 loadings")+
  scale_x_discrete(labels=c("1" = "temp", "2" = "precip", "3" = "solar", "4"="soil"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"), axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12))
dev.off()

rand_correlation_list[[1]][]
