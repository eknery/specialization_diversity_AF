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
linear_params = c()
sigm_params = c()
model_fits = list()

### setting optimization settings
control = list(method= "fftC", parscale=0.1, reltol=0.001)

### setting linear function
xr = range(trait_values) + 10*c(-se_trait, se_trait)
linear.x = make.linear.x(x0=xr[1], x1=xr[2])

for (i in 1:length(phylo_trees)){
  ### pick a phylogenetic tree
  i = 1
  one_tree = phylo_trees[[i]] 
  ### starting parameter values
  start_values = starting.point.quasse(one_tree, states=trait_values)
  ### setting quasse functions
  # constant
  const =  make.quasse(one_tree, states=trait_values, states.sd=se_trait , lambda=constant.x, mu=constant.x, sampling.f=0.9)
  #const = constrain(const, drift ~ 0)
  # linear
  linear = make.quasse(one_tree,  states=trait_values, states.sd=se_trait , lambda=linear.x, mu=constant.x, sampling.f=0.9)
  #linear = constrain(linear, drift ~ 0)
  # sigmoid
  sigm = make.quasse(one_tree, states=trait_values, states.sd=se_trait , lambda=sigmoid.x, mu=constant.x, sampling.f=0.9)
  #sigm = constrain(sigm, drift ~ 0)
  ### optimization
  ## constant
  # initial values
  init_const = c(start_values[1], start_values[2], 0.01, start_values[3])
  lower_const = c(0,0,0,0)
  names(lower_const) = names(init_const) = argnames(const)
  # finding constant mle
  mle_const = find.mle(const, x.init=init_const, lower=lower_const, control=control)
  const_params = rbind(const_params, mle_const$par)
  ## linear 
  # initial values
  init_linear = c(mle_const$par[1], lm=0.01, mle_const$par[2], mle_const$par[3:4])
  names(init_linear) = argnames(linear)
  # finding linear mle
  mle_linear = find.mle(linear, x.init=init_linear, control=control)
  linear_params = rbind(linear_params, mle_linear$par)
  ## sigmoid 
  # initial values
  init_sigm = c(mle_const$par[1], mle_const$par[1], l.xmid=mean(xr), lr=1, mle_const$par[3:4])
  names(init_sigm) = argnames(sigm)
  # finding sigmoid mle
  mle_sigm = find.mle(sigm, x.init=init_sigm, control=control)
  sigm_params= rbind(sigm_params, mle_sigm$par)
  ### summarizing model fit
  const_fit = c(lnlik= mle_const$lnLik, n_par=length(mle_const$par))
  linear_fit = c(mle_linear$lnLik, length(mle_linear$par))
  sigm_fit = c(mle_sigm$lnLik, length(mle_sigm$par))
  fit_values = rbind(const_fit, linear_fit, sigm_fit)
  model_fits[[i]] = fit_values
  ### update!
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
write.table(model_fits_df, "quasse/model_fits_df.csv", sep=",", quote=F, row.names = T)
write.table(const_params, "quasse/const_params.csv", sep=",", quote=F, row.names = F)
write.table(linear_params, "quasse/linear_params.csv", sep=",", quote=F, row.names = F)
write.table(sigm_params, "quasse/sigm_params.csv", sep=",", quote=F, row.names = F)

########################### choosing the best ###########################

model_fits_df= read.table("4_quasse/model_fits_df.csv", sep=",", h = T)

final_set= nrow(model_fits_df)-2
best_fit_per_tree =c()
for (i in seq(from=1, to=final_set, by=3)){
  one_tree = model_fits_df[i:(i+2),]
  best_fit_per_tree = rbind(best_fit_per_tree, one_tree[one_tree$aicc==min(one_tree$aicc),] )
}

 