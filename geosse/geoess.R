setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

library(ape)
library(diversitree)
library(FAmle)

### loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### define states threshold
high_ths = 0.90
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "2"
geo_states[af_percentage <= low_ths] = "1"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "0"
geo_states = as.numeric(geo_states)
names(geo_states) = spp_count_domain$species


############################### fitting models across trees  ##########################
d_params = c()
sd_params = c()
sxd_params = c()
model_fits = list()

for (i in 1:length(phylo_trees)){
  ### pick a phylogenetic tree
  one_tree = phylo_trees[[i]] 
  ### set models
  sxd = make.geosse(one_tree, states=geo_states, sampling.f=0.9)
  sd = constrain(sxd, xA ~ xB)
  d = constrain(sxd, sA ~ sB, xA ~ xB)
  ### starting values
  start_geosse = starting.point.geosse(one_tree)
  ### find mle
  # all constant
  mle_d = find.mle(d, start_geosse)
  d_params = rbind(d_params, mle_d$par)
  # varying lambda
  mle_sd = find.mle(sd, start_geosse)
  sd_params = rbind(sd_params, mle_sd$par)
  # varying lambda and mu
  mle_sxd = find.mle(sxd, start_geosse)
  sxd_params = rbind(sxd_params, mle_sxd$par)
  ### summarizing model fit
  d_fit = c(lnlik= mle_d$lnLik, n_par=length(mle_d$par))
  sd_fit = c(mle_sd$lnLik, length(mle_sd$par))
  sxd_fit = c(mle_sxd$lnLik, length(mle_sxd$par))
  fit_values = rbind(d_fit, sd_fit, sxd_fit)
  model_fits[[i]] = fit_values
  print(paste("Time:", Sys.time(), "Loop iterarion:", as.character(i) ) )
}

### arranging into dara frame
model_fits_df = data.frame()
for (i in 1:length(model_fits)){
  model_fits_df = rbind(model_fits_df, model_fits[[i]])
}

### AIC = -2(log-likelihood) + 2K
aic = -2*model_fits_df$lnlik +2*(model_fits_df$n_par)
### AiC = AIC - 2k(k+1) / (n-k-1)
aicc = aic -2*(model_fits_df$n_par+1)/(66-model_fits_df$n_par-1)
### adding
model_fits_df = data.frame(model_fits_df, aic, aicc)

### exporting
write.table(model_fits_df, "geosse/geosse_model_fits_df.csv", sep=",", quote=F, row.names = T)
write.table(d_params, "geosse/d_params.csv", sep=",", quote=F, row.names = F)
write.table(sd_params, "geosse/sd_params.csv", sep=",", quote=F, row.names = F)
write.table(sxd_params, "geosse/sxd_params.csv", sep=",", quote=F, row.names = F)

########################### choosing the best ###########################

model_fits_df= read.table("geosse/geosse_model_fits_df.csv", sep=",", h = T)

best_fit_per_tree =c()
index = seq(1, 298, by =3)
for (i in index){
  one_set = model_fits_df[i:(i+2),]
  first_lowest = one_set[one_set$aicc==min(one_set$aicc),]
  minus_first = one_set[-which(one_set$aicc==min(one_set$aicc) ),]
  second_lowest = minus_first[minus_first$aicc==min(minus_first$aicc),]
  if (second_lowest$aicc - first_lowest$aicc < 2){
    best_fit_per_tree = rbind(best_fit_per_tree, second_lowest)
  } else {
    best_fit_per_tree = rbind(best_fit_per_tree, first_lowest)
  }
}

model_name = c()
for(i in 1:nrow(best_fit_per_tree)){
  str_names = strsplit(rownames(best_fit_per_tree), '_')
  model_name = c(model_name, str_names[[i]][1])
}

best_fit_models = data.frame(model_name, best_fit_per_tree)
table(best_fit_models$model_name)

### exporting
write.table(best_fit_models, "geosse/geosse_best_fit_models.csv", sep=",", quote=F, row.names = F)


