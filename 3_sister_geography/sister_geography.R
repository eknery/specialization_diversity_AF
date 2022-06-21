setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

### packages
library(raster)
library(sp)
library(sf)
library(spatialEco)
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
