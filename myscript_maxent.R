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

############################# comparative niche analyses ##############################
library(ape)
library(phytools)
library(geiger)
library(OUwie)
library(nlme)

### loading tree
mcc=read.tree("C:/Users/eduar/Desktop/mcc_phylo.nwk")

### loading niche values
spp_niches_table=read.table("spp_niches_table.csv", header =T, sep=",",  na.strings = "NA", fill=T)
str(spp_niches_table)

### median values
median_spp_niches = aggregate(spp_niches_table[,2:5], by=list(spp_niches_table[,1]), median)
env_trait = as.matrix(median_spp_niches[,2:5])
rownames(env_trait) = median_spp_niches$Group.1

### classifying
AF=c("atlantica","baumgratziana","brunnea","budlejoides","capixaba", 
     "castaneiflora","cinerascens","discolor","dura","fasciculata",
     "formosa","hyemalis","kollmannii","kriegeriana","lymanii","mellina",
     "octopetala","penduliflora","petroniana","polyandra","racemifera",
     "robusta","ruschiana","setosociliata","shepherdii","valtheri","willdenowii")

distribution = median_spp_niches$Group.1
distribution[distribution %in% AF] = "AF-endemic"
distribution[distribution != "AF-endemic"] = "widespread"

### geographic states
states = distribution
names(states) = median_spp_niches$Group.1

### PCA
pca = prcomp(median_spp_niches[,2:5], scale=T, center = T)
pc_scores = data.frame(distribution, pca$x)
pc_trait = pca$x[,1:2]
rownames(pc_trait) = median_spp_niches$Group.1

### simmap
simmaps=make.simmap(mcc, x=states, model="ER", nsim=100)
mean_map=summary(simmaps)

# setting ancestral states
anc_states= rep(NA, length(mean_map$ace[,1]))
for (i in 1:length(mean_map$ace[,1])){
  bool = mean_map$ace[i,] == max(mean_map$ace[i,])
  high_prob_state = names(bool)[bool]
  anc_states[i] = high_prob_state
}

names(anc_states) = row.names(mean_map$ace)

# visual parameters
edge_cols = anc_states
edge_cols[edge_cols == 'AF-endemic'] = "green3"
edge_cols[edge_cols == 'widespread'] = "magenta"
node_cols = edge_cols[1:65]
tip_cols = edge_cols[66:131]
tip_cols = tip_cols[mcc$tip.label]
cols = c(tip_cols, node_cols) 
names(cols)= 1:(length(mcc$tip)+mcc$Nnode)

### phylomorphospace
phylomorphospace(mcc, env_trait[,1:2], label='off', control=list(col.node=cols))



