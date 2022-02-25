setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

### loading dataset and packages
# packages
library(raster)
library(sf)
library(sp)
library(rgeos)
library(dismo)
library(spatialEco)

library(rJava)

#loading coordinates
spp_points=read.table("spp_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)
str(spp_points)

# loading melastom_bg
bg_points=read.table("bg_points_7km.csv", header =T, sep=",",  na.strings = "NA", fill=T)
str(bg_points)

### loading raster layers
ras1 = raster("C:/Users/eduar/Desktop/rasters Neotropico 2.5/isothermality")
ras2 = raster("C:/Users/eduar/Desktop/rasters Neotropico 2.5/precipitation_seasonality")
ras3 = raster("C:/Users/eduar/Desktop/rasters Neotropico 2.5/mean_solar_radiation")
ras4 = raster("C:/Users/eduar/Desktop/rasters Neotropico 2.5/soil_pH.gri")
env_ras= stack(ras1,ras2,ras3,ras4)

alt_ras = raster("C:/Users/eduar/Desktop/rasters Neotropico 2.5/altitude")

### croping rasters
#set extent
ext=extent(-100, -33, -34, 24)
allras_crop=crop(env_ras[[1]], ext)
for (i in 2:length(env_ras@layers) ){
  ras_crop =  crop(env_ras[[i]], ext)
  allras_crop= raster::stack(allras_crop, ras_crop)
}
crop_env_ras = allras_crop


############################ MaxEnt modeling  ########################

### setting kfold validation function
maxent_kfold_evaluation = function (sp_points, bg_points, predictors, km_distance, k_value){
  performance = data.frame(matrix(NA, nrow=k_value, ncol=3))
  colnames(performance) = c('max_TSS', 'TSS_th', 'AUC')
  # preparing bg spatialpoint-format
  bg_sp = SpatialPoints(bg_points)
  crs(bg_sp) = '+proj=longlat +datum=WGS84 +no_defs'  
  # k-fold
  sp_folds = kfold(sp_points, k=k_value)
  for (i in 1:k_value){
    sp_train = sp_points[sp_folds != i,]
    sp_test = sp_points[sp_folds == i,]
    # modeling
    me = maxent(x=predictors, p=sp_train, a=bg_points, removeDuplicates=T)
    # creating pseudo-absences = all sites - sp occurrences
    circ_around <- circles(sp_train, d= km_distance, lonlat=TRUE)
    poly_around <- polygons(circ_around)
    bg_less = erase.point(bg_sp, poly_around, inside = TRUE)
    # evaluating model
    eval_me = evaluate(model=me, p=sp_test, a=bg_less, x=predictors) 
    # picking performance metrics
    th_values = eval_me@t
    tss_values = eval_me@TPR + eval_me@TNR - 1
    performance$max_TSS[i] = max(tss_values)
    performance$TSS_th[i] = th_values[tss_values == max(tss_values)][1]
    performance$AUC[i] = eval_me@auc
  }
  return(performance)
}

### looping over species names
spp_names = unique(spp_points$species)
performance_list = vector('list', length(spp_names))
names(performance_list) = spp_names

for (sp in spp_names){
  index = which(spp_names == sp)
  one_sp_points = spp_points[spp_points$species == sp, 2:3]
  try( 
       if (nrow(one_sp_points) >= 10) {
         performance_list[[index]] = maxent_kfold_evaluation(sp_points = one_sp_points, bg_points = bg_points, predictors= crop_env_ras, km_distance = 10000, k_value = 5) 
       } 
       else {
         performance_list[[index]] = maxent_kfold_evaluation(sp_points = one_sp_points, bg_points = bg_points, predictors= crop_env_ras, km_distance = 10000, k_value = nrow(one_sp_points) )
       }
  )
}

### organizing in a table
performance_table = data.frame(matrix(0,nrow=1,ncol=4))
colnames(performance_table) = c('species','max_TSS', 'TSS_th', 'AUC')

for (i in 1:length(performance_list)){
  if( is.null(performance_list[[i]]) ) {next}
  sp = names(performance_list)[i]
  species= rep(sp, length.out=nrow(performance_list[[i]]) )
  one_sp_table = data.frame(species, performance_list[[i]]) 
  performance_table = rbind(performance_table, one_sp_table)
}

performance_table = performance_table[-1,]
str(performance_table)
write.table(performance_table,"performance_table.csv", sep=",", row.names = F, quote = F)

### mean and performance table
mean_performance = aggregate(performance_table[,-1], by= list(performance_table[,1]), mean)
mean_performance
write.table(mean_performance,"mean_performance.csv", sep=",", row.names = F, quote = F)

####################### comparing niche similarity ##################

spp_names = unique(spp_points$species)
warren_table = data.frame(matrix(NA, nrow=length(spp_names), ncol=length(spp_names)))
colnames(warren_table) = spp_names
rownames(warren_table) = spp_names

for (j in 1:length(spp_names)) {
  j_name = spp_names[j]
  j_sp = spp_points[spp_points$species == j_name,2:3]
  for (i in j:length(spp_names)){
    i_name = spp_names[i]
    i_sp = spp_points[spp_points$species == i_name,2:3]
    niche_eq = nicheEquivalency(sp1=j_sp, sp2=i_sp, predictors=crop_env_ras, n=1, model=maxent)
    warren_table[i,j] = niche_eq$statistic[2]
  }
}

write.table(warren_table, "warren_table.csv", sep=',', quote= F, row.names=F)


######################### inferring environmental niche ###########################

mean_performance = read.table('mean_performance.csv', sep=',', h=T)
spp_names = unique(spp_points$species)
spp_niches_list = vector('list', length(spp_names))
names(spp_niches_list) = spp_names

for (sp_name in spp_names){
  index = which(spp_names == sp_name)
  one_sp = spp_points[spp_points$species == sp_name, 2:3]
  sp_th = mean_performance[mean_performance$Group.1 == sp_name,3]
  me = maxent(x=crop_env_ras, p=one_sp, a=bg_points, removeDuplicates=T)
  project = predict(me, crop_env_ras)
  sp_me_env = crop_env_ras[(project > sp_th)]
  spp_niches_list[[index]] = sp_me_env[sample(1:nrow(sp_me_env), 100),]
}


### organizing in a table
spp_niches_table = data.frame(matrix(0,nrow=1,ncol=5))
colnames(spp_niches_table) = c('species', 'isothermality','precip_seasonality', 'solar_radiation', 'soil_pH')

for (i in 1:length(spp_niches_list)){
  if( is.null(spp_niches_list[[1]]) ) {next}
  sp = names(spp_niches_list)[i]
  species= rep(sp, length.out=nrow(spp_niches_list[[i]]) )
  one_sp_table = data.frame(species, spp_niches_list[[i]]) 
  spp_niches_table = rbind(spp_niches_table, one_sp_table)
}

spp_niches_table = spp_niches_table[-1,]
str(spp_niches_table)

write.table(spp_niches_table, "spp_niches_table.csv", sep=',', quote= F, row.names=F)

############################# niche structure ##############################
library(tidyverse)
library(PupillometryR)
library(ggpubr)
library(readr)
library(tidyr)

### loading niche values
spp_niches_table=read.table("spp_niches_table.csv", header =T, sep=",",  na.strings = "NA", fill=T)
str(spp_niches_table)

### mean values
mean_spp_niches = aggregate(spp_niches_table[,2:5], by=list(spp_niches_table[,1]), mean)

###### PCA 
correlation = cor(mean_spp_niches[,2:5])
pca = prcomp(mean_spp_niches[,2:5], scale=T, center = T)
stdev = pca$sdev / sum(pca$sdev)
pc_number =c("1", "2", "3", "4")
stdev = data.frame(pc_number, stdev)

### randomization procedure 
rand_correlation_list = vector('list', 99)
rand_stdev_df = data.frame(matrix(0, nrow=1, ncol=4))
rand_rotation_list = vector('list', 99)

for (i in 1:99){
  isothermality = mean_spp_niches[sample(1:nrow(mean_spp_niches)), 2]
  precip_seasonality = mean_spp_niches[sample(1:nrow(mean_spp_niches)), 3]
  solar_radiation = mean_spp_niches[sample(1:nrow(mean_spp_niches)), 4]
  soil_pH = mean_spp_niches[sample(1:nrow(mean_spp_niches)), 5]
  rand_mean_niches= cbind(isothermality, precip_seasonality, solar_radiation, soil_pH)
  rand_correlation_list[[i]] = cor(rand_mean_niches)
  rand_pca = prcomp(rand_mean_niches, scale=T)
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
  geom_errorbar(aes(ymin=X5.*100, ymax=X95.*100), col="gray", size=1, width=0.25)+
  geom_point(data=rand_stdev_stats, aes(y=rand_stdev_means*100), color="gray", size = 2) +
  geom_point(data=stdev, aes(y=stdev*100), color="red", size = 2) +
  labs(x="PC axis", y="% explained variance")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12))
dev.off()
  
### collecting rotation statistics
rand_rotation_df = rand_rotation_list[[1]][,1]
for (i in 2:length(rand_rotation_list) ){
  rand_rotation_df = rbind(rand_rotation_df , rand_rotation_list[[i]][,1])
}
rand_rotation_means = apply(rand_rotation_df, 2,FUN=mean)
rand_rotation_quantiles = apply(rand_rotation_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
rand_rotation_stats= t(rbind(rand_rotation_means, rand_rotation_quantiles))
variable = c("1","2","3","4")
rand_rotation_stats = data.frame(variable, rand_rotation_stats)

# observed rotation of PC1
rotation = pca$rotation[,1]
rotation = data.frame(variable, rotation)

### ploting
tiff("PC1_loadings.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data=rand_rotation_stats, aes(x=variable, y=rand_rotation_means) ) +
  geom_errorbar(aes(ymin= X5., ymax= X95.), col="gray", size=1, width=0.25)+
  geom_point(data=rand_rotation_stats, aes(y=rand_rotation_means), color="gray", size = 2) +
  geom_point(data=rotation, aes(y=rotation), color="red", size = 2) +
  labs(x="variable", y="PC1 loadings")+
  scale_x_discrete(labels=c("1" = "therm", "2" = "precip", "3" = "solar", "4"="soil"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"), axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12))
dev.off()
