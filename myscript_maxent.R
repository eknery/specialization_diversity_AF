setwd("C:/Users/eduar/Documents/GitHub/specialization_diversification_AF")

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
# checking limits based on background
range(bg_points[,1])
range(bg_points[,2])

#set extent
ext=extent(-100, -33, -34, 24)

allras_crop=crop(env_ras[[1]], ext)

for (i in 2:length(env_ras@layers) ){
  ras_crop =  crop(env_ras[[i]], ext)
  allras_crop= raster::stack(allras_crop, ras_crop)}

crop_env_ras = allras_crop

############################ MaxEnt modeling  ########################

all_spp_names = unique(spp_points$species)
bg_sp = SpatialPoints(bg_points)
crs(bg_sp) = '+proj=longlat +datum=WGS84 +no_defs'
  
sp_name = 'discolor'
index = which(all_spp_names == sp_name)
sp_point = spp_points[spp_points$species == sp_name,2:3]
# k-fold
sp_folds = kfold(sp_point, k=5)
sp_train = sp_point[sp_folds != 1,]
sp_test = sp_point[sp_folds == 1,]
# modeling
me = maxent(x=crop_env_ras, p=sp_train, a=bg_points, removeDuplicates=T)
# creating pseudo-absences = all sites - sp occurrences
circ_around <- circles(sp_train, d=10000, lonlat=TRUE)
poly_around <- polygons(circ_around)
bg_less = erase.point(bg_sp, poly_around, inside = TRUE)
# evaluating model
eval_me = evaluate(model=me, p=sp_test, a=bg_less, x=crop_env_ras) 
th_values = eval_me@t
eval_me@auc
tss_values = eval_me@TPR + eval_me@TNR - 1
max(tss_values)
th_values[tss_values==max(tss_values)]

threshold(eval_me)



density(eval_me)
# response curves
response(me)

# projecting
project_3 = predict(me, crop_env_ras)
plot(project_3)
points(sp_train, cex=0.5)
points(bg_less, cex=0.5, col="red")

 x
### extracting suitable environmental values
th = as.numeric(threshold(eval_me)[1]) # max kappa
suitable = env_ras[project_3 > th]

hist(suitable[,1])
hist(suitable[,2])
hist(suitable[,3])
hist(suitable[,4])
