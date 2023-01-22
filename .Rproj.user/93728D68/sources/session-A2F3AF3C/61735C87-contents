
library(ape)
library(diversitree)
library(FAmle)

### loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### define states threshold
high_ths = 0.9
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "2"
geo_states[af_percentage <= low_ths] = "1"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "0"
geo_states = as.numeric(geo_states)
names(geo_states) = spp_count_domain$species


############################### fitting models across trees  ##########################
d_lin_params = c()
sd_s_lin_params = c()
sxd_sx_lin_params = c()
model_fit_list = list()

### time functions
t_functions = rep(c("linear.t", "constant.t"),c(5, 2))

for (i in 1:length(phylo_trees)){
  ### pick a phylogenetic tree
  one_tree = phylo_trees[[i]] 
  ### set models
  sxd_sx_lin = make.geosse.t(one_tree, states=geo_states, functions=t_functions, sampling.f=0.9)
  sd_s_lin = constrain(sxd_sx_lin, xA.c ~ xB.c, xA.m ~ xB.m)
  d_lin = constrain(sd_s_lin , sA.c ~ sB.c, sA.m ~ sB.m)
  ### starting values = sA.c, sA.m, sB.c, sB.m, sAB.c, sAB.m, xA, xB, dA, dB
  p = starting.point.geosse(one_tree)
  start_geosse_t = c(p[1], p[1], p[2], p[2], p[3], p[3], p[4], p[4], p[5], p[5], p[6], p[7])
  names(start_geosse_t) = argnames(sxd_sd_lin)
  ### find mle
  # lambda and mu constrained
  mle_d_lin = find.mle(d_lin, start_geosse_t)
  d_lin_params = rbind(d_lin_params, mle_d_lin$par)
  # mu constrained
  mle_sd_s_lin = find.mle(sd_s_lin, start_geosse_t)
  sd_s_lin_params = rbind(sd_s_lin_params, mle_sd_s_lin$par)
  # all free parameters
  mle_sxd_sx_lin = find.mle(sxd_sx_lin, start_geosse_t)
  sxd_sx_lin_params = rbind(sxd_sx_lin_params, mle_sxd_sx_lin$par)
  ### summarizing model fit
  d_lin_fit = c(lnlik= mle_d_lin$lnLik, n_par=length(mle_d_lin$par))
  sd_s_lin_fit = c(mle_sd_s_lin$lnLik, length(mle_sd_s_lin$par))
  sxd_sx_lin_fit = c(mle_sxd_sx_lin$lnLik, length(mle_sxd_sx_lin$par))
  fit_values = rbind(d_lin_fit, sd_s_lin_fit, sxd_sx_lin_fit)
  model_fit_list[[i]] = fit_values
  print(paste("Time:", Sys.time(), "Loop iterarion:", as.character(i) ) )
}

### arranging into dara frame
model_fit_df = data.frame()
for (i in 1:length(model_fit_list)){
  model_fit_df = rbind(model_fit_df, model_fit_list[[i]])
}

### AIC = -2(log-likelihood) + 2K
aic = (2*model_fit_df$n_par) - (2*model_fit_df$lnlik)
### AiCc = AIC - 2k(k+1) / (n-k-1)
aicc = aic + (((2*model_fit_df$n_par)*(model_fit_df$n_par+1))/(66-model_fit_df$n_par-1))
### adding
model_fit_df$aic =  aic
model_fit_df$aicc = aicc

### exporting
write.table(model_fit_df, "5_geosse_time/geosse_time_model_fit_df.csv", sep=",", quote=F, row.names = T)
write.table(d_lin_params, "5_geosse_time/d_lin_params.csv", sep=",", quote=F, row.names = F)
write.table(sd_s_lin_params, "5_geosse_time/sd_s_lin_params.csv", sep=",", quote=F, row.names = F)
write.table(sxd_sx_lin_params, "5_geosse_time/sxd_sx_lin_params.csv", sep=",", quote=F, row.names = F)

########################### choosing the best ###########################

model_fit_df= read.table("5_geosse_time/geosse_time_model_fit_df.csv", sep=",", h = T)

best_fit_per_tree =c()
index = seq(1, 298, by =3)
for (i in index){
  one_set = model_fit_df[i:(i+2),]
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
write.table(best_fit_models, "5_geosse_time/geosse_time_best_fit_models.csv", sep=",", quote=F, row.names = F)

#################### summarizing best model estimates #########################
best_fit_models = read.table("5_geosse_time/geosse_time_best_fit_models.csv", sep=",", h = T)

d_lin_params = read.table("5_geosse_time/d_lin_params.csv", sep=",",  h = T)
sd_s_lin_params = read.table("5_geosse_time/sd_s_lin_params.csv", sep=",",  h = T)
sxd_sx_lin_params = read.table("5_geosse_time/sxd_sx_lin_params.csv", sep=",", h = T)

### most common model
most_repeated = max(table(best_fit_models$model_name))
common_model = names(which(table(best_fit_models$model_name) == most_repeated) ) 
common_scenario = which(best_fit_models$model_name == common_model)
common_params = sd_s_lin_params[common_scenario,]

### descriptive statistics
apply(common_params, MARGIN = 2, FUN=mean)
apply(common_params, MARGIN = 2, FUN= function(x){ 1.96* ( sd(x)/sqrt( length(x) ) ) } )

### export best-fit model parameters
write.table(common_params, "5_geosse_time/geosse_time_best_fit_params.csv", sep=",", quote=F, row.names = F)

