setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

### packages
library(ape)
library(phytools)
library(geiger)
library(OUwie)

# loading environmental  niche
spp_niches_table=read.table("spp_niches_table.csv", header =T, sep=",",  na.strings = "NA", fill=T)
str(spp_niches_table)

### mean and se values
mean_spp_niches = aggregate(spp_niches_table[,2:5], by=list(spp_niches_table[,1]), mean)
env_traits = as.matrix(mean_spp_niches[,2:5])
rownames(env_traits) = mean_spp_niches$Group.1

se_spp_niches = aggregate(spp_niches_table[,2:5], by=list(spp_niches_table[,1]), function(x){sd(x)/sqrt(100)} )

### loading trees
mcc=read.tree("mcc_phylo.nwk")
trees=read.tree("100_rand_phylos.nwk")

### vector for classification
AF=c("atlantica","baumgratziana","brunnea","budlejoides","capixaba", 
     "castaneiflora","cinerascens","discolor","dura","fasciculata",
     "formosa","hyemalis","kollmannii","kriegeriana","lymanii","mellina",
     "octopetala","penduliflora","petroniana","polyandra","racemifera",
     "robusta","ruschiana","setosociliata","shepherdii","valtheri","willdenowii")

distribution = mean_spp_niches$Group.1
distribution[distribution %in% AF] = "AF-endemic"
distribution[distribution != "AF-endemic"] = "widespread"

### geographic states
states = distribution
names(states) = mean_spp_niches$Group.1

############################ single stochastic mapping ########################
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
phylomorphospace(mcc, env_trait[,c(1,3)], label='off', control=list(col.node=cols))

########################## fitting evolutionary models ##########################

### setting function to fit different models
fit_evo_models = function(tree, regimes, models_to_fit){
  states = regimes[,2]
  names(states) = regimes[,1]
  # reconstructing ancestral states
  simmaps=make.simmap(tree, x=states, model="ER", nsim=100)
  mean_map=summary(simmaps)
  anc_states= rep(NA, length(mean_map$ace[,1]))
  for (i in 1:length(mean_map$ace[,1])){
    bool = mean_map$ace[i,] == max(mean_map$ace[i,])
    high_prob_state = names(bool)[bool]
    anc_states[i] = high_prob_state
  }
  tree$node.label=anc_states[1:65]
  #setting fitting tables
  model_fit_table = data.frame(matrix(NA, nrow= length(models_to_fit), ncol=3))
  colnames(model_fit_table) = c("model","llik","aicc")
  #setting estimate tables
  model_estimate_list = vector("list", length(models_to_fit))
  for (i in 1:length(models_to_fit)){
    # fitting models
    fit=OUwie(tree, regimes, models_to_fit[i], mserr = "known", ub=Inf) 
    # picking fitting metrics
    model_fit_table[i,] = c(models_to_fit[i],fit$loglik,fit$AICc)
    # picking model estimates
    model_estimate = vector("list", 2)
    model_estimate[[1]] = fit$theta
    model_estimate[[2]] = fit$solution
    names(model_estimate) = c("theta","solution")
    model_estimate_list[[i]] = model_estimate
    names(model_estimate_list)[i] = models_to_fit[i]
  }
  fit_results = vector("list", 2)
  fit_results[[1]] = model_fit_table
  fit_results[[2]] = model_estimate_list
  names(fit_results) = c("fit_metrics","model_estimates")
  return(fit_results)
}

### setting function to choose best-fit model and estimates
choose_best = function (fit_results){
  if (min(fit_results$fit_metrics$aicc) < 0){
    # exclude erroneous model
    fit_results$fit_metrics = fit_results$fit_metrics[-which(fit_results$fit_metrics$aicc < 0),]
  }
  # calculate delta_aicc 
  delta_aicc = as.numeric(fit_results$fit_metrics$aicc) - as.numeric(min(fit_results$fit_metrics$aicc))
  fit_results$fit_metrics = data.frame(fit_results$fit_metrics, delta_aicc)
  # find first and second lowest delta aicc
  first_delta = fit_results$fit_metrics[fit_results$fit_metrics$delta_aicc == min(fit_results$fit_metrics$delta_aicc),]
  minus_first = fit_results$fit_metrics[-which(fit_results$fit_metrics$delta_aicc == min(fit_results$fit_metrics$delta_aicc)),]
  second_delta =  minus_first[minus_first$delta_aicc == min(minus_first$delta_aicc),]
  # index to check complexity
  model_names = fit_results$fit_metrics[,1]
  first_index = which(model_names == first_delta$model)
  second_index = which(model_names == second_delta$model)
  # compare delta aicc and pick best-fit model
  if(first_index > second_index){
    if (second_delta$delta_aicc - first_delta$delta_aicc > 2){
      best_fit_model = first_delta
    } else {
      best_fit_model = second_delta
    }
  } else {
    best_fit_model = first_delta
  }
  # finding best-model estimates
  best_fit_estimates = fit_results$model_estimates[names(fit_results$model_estimates) == best_fit_model$model]
  best_fit_estimates = best_fit_estimates[[1]]
  # pick best-fit estimates
  theta_se = c()
  for(j in 1:nrow(best_fit_estimates$theta) ){
    thse = c(best_fit_estimates$theta[j,])
    names(thse)[1] = paste("theta", as.character(j), sep="_")
    names(thse)[2] = paste("se", as.character(j), sep="_")
    theta_se = c(theta_se, thse)
  }
  alp_sig = c()
  for(j in 1:ncol(best_fit_estimates$solution) ){
    alsi=best_fit_estimates$solution[,j]
    names(alsi)[1] = paste("alpha", as.character(j), sep="_")
    names(alsi)[2] = paste("sigma_sq", as.character(j), sep="_")
    alp_sig = c(alp_sig, alsi)
  }
  best_estimate = c(theta_se,alp_sig)
  # packing things up
  best_fit_metric = best_fit_model[1,]
  best_list = list(best_fit_metric[-5], best_estimate)
  names(best_list) = c("best_fit", "best_estimates")
  return(best_list)
}

### fitting and choosing models across trees
# setting regimes & models
regimes = data.frame(mean_spp_niches$Group.1, distribution, mean_spp_niches$solar_radiation, se_spp_niches$solar_radiation)
colnames(regimes) = c("species","state", "trait", "se") 
models=c("BM1","BMS","OU1","OUM","OUMV")

# looping across trees
all_best_models = data.frame(matrix(NA, nrow=length(trees), ncol=4))
colnames(all_best_models) = c(c("model","llik","aicc","delat_aicc"))
all_best_estimates = vector("list" , length(trees))

for (i in 1:length(trees)){
  one_tree = trees[[i]]
  all_fit = fit_evo_models(tree= one_tree, regimes= regimes, models_to_fit= models)
  best_choice = choose_best(all_fit)
  all_best_models[i,] = best_choice$best_fit
  all_best_estimates[[i]] = best_choice$best_estimates
}

### organizing into a dataframe
all_best_estimates_df = t(data.frame(all_best_estimates[[1]]))
for(i in 2:length(all_best_estimates)){
  all_best_estimates_df = rbind(all_best_estimates_df,  t(data.frame(all_best_estimates[[i]])) )
}
all_best_estimates_df = data.frame(all_best_models$model, all_best_estimates_df)
colnames(all_best_estimates_df)[1] = "model"

#exporting
write.table(all_best_models, "solar_best_models.csv", sep=',', quote=F, row.names = F)
write.table(all_best_estimates_df, "solar_best_estimates.csv", sep=',', quote=F, row.names = F)

############################ estimating neutral evolution ######################

regimes = data.frame(mean_spp_niches$Group.1, distribution, mean_spp_niches$precip_seasonality, se_spp_niches$precip_seasonality)
colnames(regimes) = c("species","state", "trait", "se") 

### fitting BM model
bm_fits=vector(mode = "list", length = length(trees))

for (i in 1:length(trees)){
  one_tree=trees[[i]]
  simmaps=make.simmap(one_tree, x=states, model="ER", nsim=100)
  mean_map=summary(simmaps)
  anc_states= rep(NA, length(mean_map$ace[,1]))
  for (j in 1:length(mean_map$ace[,1])){
    bool = mean_map$ace[j,] == max(mean_map$ace[j,])
    high_prob_state = names(bool)[bool]
    anc_states[j] = high_prob_state
  }
  one_tree$node.label=anc_states[1:65]
  fit=OUwie(one_tree, regimes, model="BM1",  mserr = "known", ub=Inf)
  bm_fits[[i]]=fit
}

### collecting parameter estimates
bm_estimates=as.data.frame(matrix(NA,nrow=length(bm_fits),ncol=5))
colnames(bm_estimates)=c("model","theta","se","alpha","sigma_q")

for (i in 1:length(bm_fits)){
  bm_estimates[i,1]=bm_fits[[i]]$model
  bm_estimates[i,2:3]=bm_fits[[i]]$theta[1,] # optima & se
  bm_estimates[i,4:5]=bm_fits[[i]]$solution[,1] # alpha & sigma_q
}

write.table(bm_estimates,"precip_bm_estimates.csv", sep=",", quote=F, row.names=F)

### simulating evolution under BM
bm_simulations=c("AF-enfdemic","widespread")

for(i in 1:length(trees)){
  one_tree=trees[[i]]
  for (j in 1:100){
    sim_hypervolumes = fastBM(one_tree, a=sample(regimes[,3],1), sig2=bm_estimates$sigma_q[i], bounds=c(0.0001,Inf), internal=F)
    af_mean= mean(sim_hypervolumes[which(names(sim_hypervolumes)%in% AF)])
    ws_mean= mean(sim_hypervolumes[-which(names(sim_hypervolumes)%in% AF)])
    bm_simulations = rbind(bm_simulations, c(af_mean, ws_mean))
  }
}

bm_simulations = data.frame(bm_simulations[-1,])
colnames(bm_simulations)<-c("AF-endemic","widespread")

write.table(bm_simulations,"precip_bm_simulations.csv", sep=",", quote=F, row.names=F)

