setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

library(ape)
library(diversitree)
library(FAmle)

### loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

### loading trait dataset
spp_hvolumes = read.table("1_hypervolume_inference/spp_hvolumes.csv", h=T, sep=",")

### setting trait vector
hvolumes = spp_hvolumes$hvolume
names(hvolumes) = spp_hvolumes$specie
trait_values= hvolumes

### sampling error
se_trait = sd(trait_values)/length(trait_values)

############################### fitting models across trees  ##########################
const_params = c()
l_lin_params = c()
lm_lin_params = c()
model_fit_list = list()

### setting optimization settings
control = list(parscale=0.1, reltol=0.001)

for (i in 1:length(phylo_trees)){
  ### pick a phylogenetic tree
  one_tree = phylo_trees[[i]] 
  ### starting parameter values and linear function
  start_quasse = starting.point.quasse(one_tree, states=trait_values)
  xr = c(range(trait_values) + c(-1,1)*start_quasse["diffusion"])
  linear.x = make.linear.x(x0=xr[1], x1=xr[2])
  ### setting quasse functions
  # constant
  const =  constrain(make.quasse(one_tree, states=trait_values, states.sd=se_trait , lambda=constant.x, mu=constant.x, sampling.f=0.9), drift~0)
  # linear lambda
  l_lin = constrain(make.quasse(one_tree,  states=trait_values, states.sd=se_trait , lambda=linear.x, mu=constant.x, sampling.f=0.9), drift~0)
  # linear lambda and mu
  lm_lin = constrain(make.quasse(one_tree, states=trait_values, states.sd=se_trait , lambda=linear.x, mu=linear.x, sampling.f=0.9), drift~0)
  ### optimization
  ## constant
  # initial values
  init_const = c(start_quasse[1], start_quasse[2], start_quasse[3])
  lower_const = c(0,0,0)
  names(lower_const) = names(init_const) = argnames(const)
  # finding constant mle
  mle_const = find.mle(const, x.init=init_const, lower=lower_const, control=control)
  const_params = rbind(const_params, mle_const$par)
  ## linear lambda 
  # initial values
  init_l_lin = c(mle_const$par[1], 0.01, mle_const$par[2], mle_const$par[3])
  names(init_l_lin) = argnames(l_lin)
  # finding mle 
  mle_l_lin = find.mle(l_lin, x.init=init_l_lin, control=control)
  l_lin_params = rbind(l_lin_params, mle_l_lin$par)
  ## linear lambda and mu
  #initial values
  init_lm_lin = c(mle_const$par[1], 0.01, mle_const$par[2], 0.01*mle_const$par[2], mle_const$par[3])
  names(init_lm_lin) = argnames(lm_lin)
  # finding mle
  mle_lm_lin = find.mle(lm_lin, x.init=init_lm_lin, control=control)
  lm_lin_params= rbind(lm_lin_params, mle_lm_lin$par)
  ### summarizing model fit
  const_fit = c(lnlik= mle_const$lnLik, n_par=length(mle_const$par))
  l_lin_fit = c(mle_l_lin$lnLik, length(mle_l_lin$par))
  lm_lin_fit = c(mle_lm_lin$lnLik, length(mle_lm_lin$par))
  fit_values = rbind(const_fit, l_lin_fit, lm_lin_fit)
  model_fit_list[[i]] = fit_values
  ### update!
  print(paste("Time:", Sys.time(), "Loop iterarion:", as.character(i) ) )
}

### arranging into dara frame
model_fit_df = data.frame()
for (i in 1:length(model_fits)){
  model_fit_df = rbind(model_fit_df, model_fit_list[[i]])
}

### AIC = -2(log-likelihood) + 2K
aic = -2*model_fit_df$lnlik +2*(model_fit_df$n_par)
### AiCc = AIC - 2k(k+1) / (n-k-1)
aicc = aic -2*(model_fit_df$n_par+1)/(66-model_fit_df$n_par-1)
### adding
model_fit_df = data.frame(model_fit_df, aic, aicc)

### exporting
write.table(model_fit_df, "quasse/quasse_model_fit_df.csv", sep=",", quote=F, row.names = T)
write.table(const_params, "quasse/const_params.csv", sep=",", quote=F, row.names = F)
write.table(l_lin_params, "quasse/l_lin_params.csv", sep=",", quote=F, row.names = F)
write.table(lm_lin_params, "quasse/lm_lin_params.csv", sep=",", quote=F, row.names = F)

########################### choosing the best ###########################

model_fit_df= read.table("quasse/quasse_model_fit_df.csv", sep=",", h = T)

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
write.table(best_fit_models, "quasse/quasse_best_fit_models.csv", sep=",", quote=F, row.names = F)

#################### summarizing best model estimates #########################
best_fit_models = read.table("quasse/quasse_best_fit_models.csv", sep=",", h = T)

const_params = read.table("quasse/const_params.csv", sep=",",  h = T)
l_lin_params = read.table("quasse/l_lin_params.csv", sep=",",  h = T)
lm_lin_params = read.table("quasse/lm_lin_params.csv", sep=",", h = T)


lm_lin_params

### most common model
most_repeated = max(table(best_fit_models$model_name))
common_model = names(which(table(best_fit_models$model_name) == most_repeated) ) 
common_scenario = which(best_fit_models$model_name == common_model)
common_params = lm_lin_params[common_scenario,]

### descriptive statistics
apply(common_params, MARGIN = 2, FUN=mean)
apply(common_params, MARGIN = 2, FUN=sd)

################################## plotting #############################
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

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### define geographic states 
high_ths = 0.9
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "AF"
geo_states[af_percentage <= low_ths] = "other"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "AFother"
names(geo_states) = spp_count_domain$species

### summary hypervolume per geographic state
aggregate(hvolumes, by=list(geo_states), mean)

### common values
x1 = rep(0, nrow(common_params))
xend = rep(max(hvolumes), nrow(common_params))
y_lim=c(0,1)

# axis names
x_axis_name = "hypervolume size"
y_axis_name = "lambda"

### line estimates
intercepts = common_params$l.c 
angulars = common_params$l.m
lines = data.frame(x1, xend, y1= intercepts, yend= xend*angulars + intercepts)

median_intercept = median(intercepts)
median_angular = median(angulars)
median_line = data.frame(x1, xend, y1= median_intercept, yend= xend*median_angular + median_intercept)

lm_lambda = ggplot() +
  geom_segment(data = lines, color="gray", size=1.2, alpha=0.1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = median_line, color="black", size=1, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  ylim(y_lim)+
  xlab(x_axis_name) + ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=9,face="bold"),axis.text=element_text(size=6),legend.position = "none")

tiff("quasse/quasse_lambda.tiff", units="in", width=2, height=1.75, res=600)
  lm_lambda
dev.off()

 