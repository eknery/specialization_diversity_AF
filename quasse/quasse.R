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
hist(sqrt(hvolumes), breaks=10)
shapiro.test(sqrt(hvolumes))
trait_values= sqrt(hvolumes)


### sampling error
se_trait = sd(trait_values)/length(trait_values)

############################### fitting models across trees  ##########################
const_params = c()
linear_params = c()
sigm_params = c()
model_fits = list()

### stopped at 24, start from 25 !!!!

for (i in 1:length(phylo_trees)){
  ### pick a phylogeentic tree
  one_tree = phylo_trees[[1]] # !!!!!!!!
  ### starting parameter values
  start_values = starting.point.quasse(one_tree, states=trait_values)
  ### setting linear function
  xr = c( range(trait_values)[1]/start_values["diffusion"],  range(trait_values)[2]+start_values["diffusion"] )
  linear.x = make.linear.x(x0=xr[1], x1=xr[2])
  ### setting quasse functions
  # constant
  quasse_const =  make.quasse(one_tree, states=trait_values, states.sd=se_trait , lambda=constant.x, mu=constant.x, sampling.f=0.9)
  quasse_const = constrain(quasse_const, drift ~ 0)
  # linear
  quasse_linear = make.quasse(one_tree,  states=trait_values, states.sd=se_trait , lambda=linear.x, mu=constant.x, sampling.f=0.9)
  quasse_linear = constrain(quasse_linear, drift ~ 0)
  # sigmoid
  quasse_sigm = make.quasse(one_tree, states=trait_values, states.sd=se_trait , lambda=sigmoid.x, mu=constant.x, sampling.f=0.9)
  quasse_sigm = constrain(quasse_sigm, drift ~ 0)
  ### setting and optimizing functions
  ## optimization control
  control = list(parscale=.1, reltol=0.001)
  ## constant
  # initial values
  init_const = c(start_values[1], start_values[2], start_values[3])
  lower_const = c(0,0,0)
  names(lower_const) = names(init_const) = argnames(quasse_const)
  # finding constant mle
  const_mle = find.mle(quasse_const, x.init=init_const, lower=lower_const, control=control)
  const_params= rbind(const_params, const_mle$par)
  ## linear 
  # initial values
  init_linear = c(const_mle$par[1], lm=0, const_mle$par[2:3])
  names(init_linear) = argnames(quasse_linear)
  # finding linear mle
  linear_mle = find.mle(quasse_linear, x.init=init_linear, control=control)
  linear_params= rbind(linear_params, linear_mle$par)
  ## sigmoid 
  # initial values
  init_sigm =  c(const_mle$par[1], const_mle$par[1], l.xmid=mean(xr), lr=1, const_mle$par[2:3])
  names(init_sigm) = argnames(quasse_sigm)
  # finding sigmoid mle
  sigm_mle = find.mle(quasse_sigm, x.init=init_sigm, control=control)
  sigm_params= rbind(sigm_params, sigm_mle$par)
  ### summarizing model fit
  const= c(lnlik= const_mle$lnLik, n_par=length(const_mle$par))
  linear= c(linear_mle$lnLik, length(linear_mle$par))
  sigm= c(sigm_mle$lnLik, length(sigm_mle$par))
  fit_values = rbind(const, linear, sigm)
  model_fits[[i]] = fit_values
  ### update!
  print(paste("Time:", Sys.time(), "Loop iterarion:", as.character(i) ) )
}


### arranging into dara frame
model_fits_df = data.frame()
for (i in 1:length(model_fits)){
  model_fits_df = rbind(model_fits_df, model_fits[[i]])
}

write.table(model_fits_df, "quasse/1_24_model_fits_df.csv", sep=",", quote=F, row.names = T)

write.table(const_params, "quasse/const_params.csv", sep=",", quote=F, row.names = F)
write.table(linear_params, "quasse/linear_params.csv", sep=",", quote=F, row.names = F)
write.table(sigm_params, "quasse/sigm_params.csv", sep=",", quote=F, row.names = F)

########################### choosing the best ###########################

model_fits_df= read.table("4_quasse/1_24_model_fits_df.csv", sep=",", h = T)

### AIC = -2(log-likelihood) + 2K
aic = -2*model_fits_df$lnlik +2*(model_fits_df$n_par)
### AIC - 2k(k+1) / (n-k-1)
aicc = aic -2*(model_fits_df$n_par+1)/(66-model_fits_df$n_par-1)
### adding
model_fits_df = data.frame(model_fits_df, aic, aicc)

final_set= nrow(model_fits_df)-2

best_fit_per_tree =c()
for (i in seq(from=1, to=final_set, by=3)){
  one_tree = model_fits_df[i:(i+2),]
  best_fit_per_tree = rbind(best_fit_per_tree, one_tree[one_tree$aicc==min(one_tree$aicc),] )
}

 