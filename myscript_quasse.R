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
names(hvolumes) = spp_hvolumes$species

### testing trait normality
hist(sqrt(hvolumes))
shapiro.test(sqrt(hvolumes))

### sampling error
sq_hvolume = sqrt(hvolumes)
se_hvolume = sd(sq_hvolume)/length(sq_hvolume)

### expected effect of independent variable
exp_effect = max(sq_hvolume) - min(sq_hvolume) / min(sq_hvolume)

### setting linear function

xr = c( range(hvolumes)[1]/start_values["diffusion"],  range(hvolumes)[2]*start_values["diffusion"] )
linear.x = make.linear.x(x0=xr[1], x1=xr[2])

############################### fitting models across trees  ##########################
const_params = c()
linear_params = c()
sigm_params = c()
model_fits = list()

for (i in 1:length(phylo_trees)){
  ### pick a phylogeentic tree
  one_tree = phylo_trees[[i]]
  ### starting parameter values
  start_values = starting.point.quasse(one_tree, states=hvolumes)
  ### setting quasse functions
  # constant
  quasse_const =  make.quasse(one_tree, states=hvolumes, states.sd=se_hvolume , lambda=constant.x, mu=constant.x, sampling.f=0.9)
  quasse_const = constrain(quasse_const, drift ~ 0)
  # linear
  quasse_linear = make.quasse(one_tree,  states=hvolumes, states.sd=se_hvolume , lambda=linear.x, mu=constant.x, sampling.f=0.9)
  quasse_linear = constrain(quasse_linear, drift ~ 0)
  # sigmoid
  quasse_sigm = make.quasse(one_tree, states=hvolumes, states.sd=se_hvolume , lambda=sigmoid.x, mu=constant.x, sampling.f=0.9)
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
  sigm= c(sigm_mle$lnLik, length(const_mle$par))
  fit_values = rbind(const, linear, sigm)
  model_fits[[1]] = fit_values
  ### update!
  print(paste("Loop iterarion:", as.character(i) ) )
}

######### visualizing #################################

y0 = sigm_mle$par[1]
y1 = sigm_mle$par[2]
xmid = sigm_mle$par[3]
r = sigm_mle$par[4]
plot(hvolumes, y0 + (y1-y0)/(1 + exp(r*(xmid-hvolumes)) ) )


