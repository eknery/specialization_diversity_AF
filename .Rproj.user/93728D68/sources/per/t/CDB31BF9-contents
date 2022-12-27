library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
library(phytools)
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

### loading mcc phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")

### number of alternative phylpogenetic trees
n_phylo = length(list.files("0_data/100_trees"))

### loading species' altitude
spp_altitude = read.table("0_data/spp_altitude.csv", sep=',', h=T)

### loading species' hypervolumes
spp_hvolumes = read.table("1_hypervolume_inference/spp_hvolumes.csv", h=T, sep=",")

# named vector of hypervolumes
hvolumes = spp_hvolumes$hvolume
names(hvolumes) = spp_hvolumes$species

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

# define geographic states 
high_ths = 0.9
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "AF"
geo_states[af_percentage <= low_ths] = "other"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "AFother"
names(geo_states) = spp_count_domain$species

# my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")
names(mycols) = c("AF", "AFother", "other")

################################## DEC reconstruction ################################

### create directory for DEC results
# check if dir exists
dir_check = dir.exists(paths="7_graphs/DEC_ancestral_reconstructions")
# create dir if not created yet
if (dir_check == FALSE){
  dir.create(path= "7_graphs/DEC_ancestral_reconstructions", showWarnings = , recursive = FALSE, mode = "0777")
}

### setting DEC model
# reading range data
geog_fn = ("0_data/spp_distribution_af.data")
moref(geog_fn)
# converting phylip format to tipranges
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geog_fn)
tipranges
# setting maximum number of areas occupied for reconstructions
max_range_size = max(rowSums(dfnums_to_numeric(tipranges@df)))
# Initialize DEC model
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
# location of the geography text file
BioGeoBEARS_run_object$geogfn = geog_fn
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size
# Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$min_branchlength = 0.001    
# set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$include_null_range = FALSE    
# computing options
BioGeoBEARS_run_object$num_cores_to_use = 1
# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

### fitting DEC over mcc tree
# phylogeny tree location
trfn = paste("0_data/mcc_phylo.nwk")
tr = read.tree(trfn)
# n tips and nodes
n_tips = Ntip(tr)
n_nodes = tr$Nnode
# inputting tree into DEC
BioGeoBEARS_run_object$trfn = trfn
# fitting DEC
res_DEC = bears_optim_run(BioGeoBEARS_run_object)
# node marginal ML
relprobs_matrix = res_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
write.table(relprobs_matrix , "7_graphs/relprobs_matrix.csv", sep=",", row.names=F, quote=F)

### plotting DEC over mcc treee
## loading relative probabilities at notes
relprobs_matrix = read.table("7_graphs/relprobs_matrix.csv", sep=",", h=T)

## setting states
# tip states probs
tip_states_probs = relprobs_matrix [1:n_tips, ]
# ancestral state probs
inner_node_probs = relprobs_matrix [(1+n_tips):(n_tips+n_nodes),]
# state colors
state_cols=c( "#1E88E5",  "#D81B60", "#FFC107")
names(state_cols)=c("AF",   "AFother", "other")

## plot
tiff("7_graphs/dec_mcc_ranges.tiff", units="cm", width=7, height=18, res=600)
  plotTree(tree=tr,fsize=0.75, ftype="off")
  tiplabels(pie=tip_states_probs, piecol=state_cols, cex=0.75)
  nodelabels(node=(1+n_tips):(n_tips+n_nodes), pie= inner_node_probs, piecol=state_cols, cex=1.5)
  axisPhylo(pos=c(0.5), font=3, cex.axis=0.5)
dev.off()

### fitting DEC over trees 
for (i in 1:n_phylo){ 
  # phylogeny tree location
  trfn = paste("0_data/100_trees/tree_", as.character(i), sep="")
  tr = read.tree(trfn)
  # n tips and nodes
  n_tips = Ntip(tr)
  n_nodes = tr$Nnode
  # inputting tree into DEC
  BioGeoBEARS_run_object$trfn = trfn
  # fitting DEC
  res_DEC = bears_optim_run(BioGeoBEARS_run_object)
  # node marginal ML
  relprobs_matrix = res_DEC$ML_marginal_prob_each_state_at_branch_top_AT_node
  # node states
  state_labels=c("AF", "other", "AFother")
  node_states = get_ML_states_from_relprobs(relprobs=relprobs_matrix, statenames=state_labels, returnwhat = "states", if_ties = "takefirst")
  # getting only ancestral nodes
  anc_node = (n_tips+1):(n_tips+n_nodes)
  state = node_states[anc_node]
  anc_node_states = data.frame(anc_node, state)
  write.table(anc_node_states , paste("7_graphs/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_"), sep=",", row.names=F, quote=F)
}

################################## LTT plot ####################################

### joining reconstruction across trees 
anc_data = c()
for (i in 1:n_phylo){ 
  # phylogeny tree location
  trfn = paste("0_data/100_trees/tree_", as.character(i), sep="")
  tr = read.tree(trfn)
  # n tips and nodes
  n_tips = Ntip(tr)
  n_nodes = tr$Nnode
  # node ages
  node_ages = round(node.depth.edgelength(tr),5)
  present = round(max(node_ages), 5)
  node_ages = node_ages - present
  # ancestral node ages
  anc_node_ages = node_ages[(n_tips+1):(n_tips+n_nodes)]
  # ancestral biogeographic state
  dec_fn = paste("7_graphs/DEC_ancestral_reconstructions/anc_node_states",as.character(i), sep="_")
  anc_node_states = read.table(dec_fn, sep=",", h=T)
  # join into dataframe
  one_anc_datum = data.frame(anc_node_states, anc_node_ages)
  # update!
  anc_data = rbind(anc_data, one_anc_datum)
}

### dividing into time intervals
# time boundaries
old_age = round(min(anc_data$anc_node_ages),2)
new_age = 0
# intervals
intervals = anc_data$anc_node_ages
breaks = round(seq(new_age, old_age, by= (old_age - new_age)/8),2)
for (i in 1:length(breaks)){
  intervals[which(anc_data$anc_node_ages >= breaks[i+1] & anc_data$anc_node_ages < breaks[i])] = round((breaks[i] + breaks[i+1])/2, 2)
}
anc_data$intervals = intervals

### cumulative lineage count per group across tree
initial_indexes = seq(1 , nrow(anc_data), by=65)
all_tree_sums = c()
for (i in initial_indexes){
  one_tree_data = anc_data[i:(i+64),]
  split_by_state = split(one_tree_data, f= one_tree_data$state)
  all_sums_df = c()
  for (ii in 1:length(split_by_state)){
    count = table(split_by_state[[ii]]$intervals)  #
    cum_sum = cumsum(count)
    state = rep(names(split_by_state)[[ii]], length(count))
    inter_age = names(count)
    one_state_df= data.frame(state, cum_sum, inter_age)
    all_sums_df = rbind(all_sums_df, one_state_df)
  }
  all_tree_sums = rbind(all_tree_sums, all_sums_df)
}

### summarize cumulative counts across trees
# split by state
split_by_state = split(all_tree_sums, f=all_tree_sums$state)
# empty object to receive final results
ltt_df = c()
for (i in 1:length(split_by_state)){
  state_data = split_by_state[[i]]
  state_name = names(split_by_state)[i]
  n_estimates = aggregate(state_data$cum_sum, by=list(state_data$inter_age), length)
  mean_count = aggregate(state_data$cum_sum, by=list(state_data$inter_age), mean)
  sd_count = aggregate(state_data$cum_sum, by=list(state_data$inter_age), sd)
  state = rep(state_name, nrow(n_estimates))
  one_state_df = data.frame(state,n_estimates,mean_count$x, sd_count$x)
  ltt_df  = rbind(ltt_df , one_state_df)
}

# naming some columns
colnames(ltt_df)[2:5] = c("inter_age", "n", "central", "dispersion")
ltt_df$inter_age = as.numeric(ltt_df$inter_age) 

# creating confidence intervals
ci = 1.96 * (ltt_df$dispersion / sqrt(ltt_df$n))
ltt_df$ci = ci

### plotting by interval
# text size
axis_title_size = 10
x_text_size = 8

ltt_plot = ggplot(data= ltt_df, aes(x=inter_age, y=central, group= state, color=state) ) +
  geom_point(size = 1, alpha = 1) +
  geom_line(size=1)+
  geom_errorbar(size=1, width=0, aes(ymin=central-ci, ymax=central+ci))+
  geom_vline(xintercept = -2.58,linetype="dotted", colour= "black", size=0.5)+
  scale_colour_manual(values=mycols)+
  xlab("time before present (m.y.a.)")+ ylab("n lineages")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),axis.text.y = element_text(angle = 90),legend.position = "none")

tiff("7_graphs/ltt_plot.tiff", units="cm", width=7, height=6, res=600)
  ltt_plot
dev.off()


####################### Describing altitude and niche breadth #############

### loading species altitude and niche breadth
spp_altitude = data.frame(geo_states, spp_altitude)
spp_hvolumes = data.frame(geo_states, spp_hvolumes)

axis_title_size = 10
x_text_size = 8

alt_plot = ggplot(data= spp_altitude, aes(x=geo_states, y=altitude, fill= geo_states)) +
  geom_point(aes(color=geo_states),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.50)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("species' elevation (m)")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),axis.text.y = element_text(angle = 90),legend.position = "none")

hv_plot = ggplot(data= spp_hvolumes, aes(x=geo_states, y=hvolume, fill= geo_states)) +
  geom_point(aes(color=geo_states),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.50)+
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("hypervolume size")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),legend.position = "none")

tiff("7_graphs/species_data.tiff", units="cm", width=7, height=12, res=600)
  ggarrange(alt_plot,hv_plot, nrow=2,ncol=1)
dev.off()

################################## RO and NO intercepts #########################

### loading geographic and niche sister data
sister_no_metrics = read.table("2_sister_hypervolume/sister_no_metrics.csv", sep=',', h=T)
sister_ro_metrics = read.table("3_sister_geography/sister_ro_metrics.csv", sep=',', h=T)

# text size
axis_title_size = 12
x_text_size = 8

### sister ro metrics
ro_plot = ggplot(data= sister_ro_metrics, aes(x=state, y=intercept_ro, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.50) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  xlab("geographic distribution")+ ylab("RO intercept")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),legend.position = "none")

### sister no metrics
no_plot = ggplot(data= sister_no_metrics, aes(x=state, y=intercept_no, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.50) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  xlab("geographic distribution")+ ylab("NO intercept")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=x_text_size),legend.position = "none")

tiff("7_graphs/intercepts_data.tiff", units="cm", width=14, height=6.5, res=600)
  ggarrange(ro_plot,no_plot, nrow=1,ncol=2)
dev.off()

######################## GeoSSE estimates ###################################

### loading fit values and model parameters
best_fit_models = read.table("4_geosse/geosse_best_fit_models.csv", sep=",", h = T)
sd_params = read.table("4_geosse/sd_params.csv", sep=",",  h = T)

### most common model
most_repeated = max(table(best_fit_models$model_name))
common_model = names(which(table(best_fit_models$model_name) == most_repeated) ) 
common_scenario = which(best_fit_models$model_name == common_model)
common_params = sd_params[common_scenario,]

### organizing estimates
state = c(rep("other", length(common_params$sA)),rep("AF", length(common_params$sB)))
s_estimates = c(common_params$sA, common_params$sB)
s_df = data.frame(state, s_estimates)

### summarizing estimates
ci_calculator = function(x){ 1.96* ( sd(x)/sqrt( length(x) ) ) }
mean_s = aggregate(s_df$s_estimates, by= list(s_df$state), mean)
ci_s = aggregate(s_df$s_estimates, by= list(s_df$state), ci_calculator)
mean_s_df = data.frame(mean_s, ci_s$x)
colnames(mean_s_df) = c("state", "mean_s", "ci_s")

### plotting
# axis title size
axis_title_size = 10
# axis names
x_axis_name = "geographic distribution"
y_axis_name = "lambda"

geo_plot= ggplot(data=mean_s_df, aes(color=state)) +
  geom_point(data=mean_s_df,mapping=aes(x=state, y=mean_s)) +
  geom_errorbar(data=mean_s_df, size=1, width=0.1,mapping=aes(x=state, ymin=mean_s-ci_s, ymax=mean_s+ci_s))+
  scale_colour_manual(values=mycols[-2])+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "other" = "outside AF"))+
  xlab(x_axis_name)+ ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size ,face="bold"),axis.text=element_text(size=6),legend.position = "none")

tiff("7_graphs/geosse_estimates.tiff", units="cm", width=7, height=6, res=600)
  geo_plot
dev.off()

####################### time-GeoSSE estimates #################################

### loading fit values and model parameters
best_fit_models = read.table("5_geosse_time/geosse_time_best_fit_models.csv", sep=",", h = T)
sd_s_lin_params = read.table("5_geosse_time/sd_s_lin_params.csv", sep=",",  h = T)

### most common model
most_repeated = max(table(best_fit_models$model_name))
frequent_model = names(which(table(best_fit_models$model_name) == most_repeated) ) 
frequent_model_n = which(best_fit_models$model_name == frequent_model)
frequent_model_params = sd_s_lin_params[frequent_model_n,]

### common values
x1 = rep(0, nrow(frequent_model_params))
xend = rep(3.5, nrow(frequent_model_params))

### parameters per group

## sA line estimates = other
intercepts_a = frequent_model_params$sA.c 
angulars_a = frequent_model_params$sA.m
lines_a = data.frame(x1, xend, y1= intercepts_a, yend = xend*angulars_a + intercepts_a)
# mean line
mean_line_a = apply(lines_a, MARGIN = 2, FUN= mean)
# confidence interval lines
ci_line_a = apply(lines_a, MARGIN = 2, FUN= sd) /  sqrt(nrow(lines_a))*1.96
ci_line_a_up = mean_line_a + ci_line_a
ci_line_a_dw =  mean_line_a - ci_line_a
# dataframe mean line est
lines_a_df = data.frame(rbind(mean_line_a, ci_line_a_up, ci_line_a_dw))

## sB line estimates = AF
intercepts_b = frequent_model_params$sB.c 
angulars_b = frequent_model_params$sB.m
lines_b = data.frame(x1, xend, y1= intercepts_b, yend = xend*angulars_b + intercepts_b)
# mean line
mean_line_b = apply(lines_b, MARGIN = 2, FUN= mean)
# confidence interval lines
ci_line_b = apply(lines_b, MARGIN = 2, FUN= sd) /  sqrt(nrow(lines_b))*1.96
ci_line_b_up = mean_line_b + ci_line_b
ci_line_b_dw =  mean_line_b - ci_line_b
# dataframe mean line est
lines_b_df = data.frame(rbind(mean_line_b, ci_line_b_up, ci_line_b_dw))

## sAB line estimates = AF and other
intercepts_ab = frequent_model_params$sAB.c 
angulars_ab = frequent_model_params$sAB.m
lines_ab = data.frame(x1, xend, y1= intercepts_ab, yend = xend*angulars_ab + intercepts_ab)
# mean line
mean_line_ab = apply(lines_ab, MARGIN = 2, FUN= mean)
# confidence interval lines
ci_line_ab = apply(lines_ab, MARGIN = 2, FUN= sd) /  sqrt(nrow(lines_ab))*1.96
ci_line_ab_up = mean_line_ab + ci_line_ab
ci_line_ab_dw =  mean_line_ab - ci_line_ab
# dataframe mean line est
lines_ab_df = data.frame(rbind(mean_line_ab, ci_line_ab_up, ci_line_ab_dw))

# make time negative
lines_a_df$xend = lines_b_df$xend = lines_ab_df$xend  = rep(-3.5,length(lines_a_df$xend))

### plotting
# axis title size
axis_title_size = 10
# axis names
x_axis_name = "time before present (m.y.a.)"
y_axis_name = "lambda"

geo_time_plot = ggplot() +
  geom_segment(data = lines_a_df[1,], color=mycols[3], size=0.75, alpha=1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = lines_a_df[2,], color=mycols[3], size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_segment(data = lines_a_df[3,], color=mycols[3], size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_segment(data = lines_b_df[1,], color=mycols[1], size=0.75, alpha=1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = lines_b_df[2,], color=mycols[1], size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_segment(data = lines_b_df[3,], color=mycols[1], size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_vline(aes(xintercept = c(-0.05,-2.58)),linetype="dotted",colour="gray",size=0.75)+
  xlab(x_axis_name)+ ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size ,face="bold"),axis.text=element_text(size=6),legend.position = "none")

tiff("7_graphs/geosse_time_estimates.tiff", units="cm", width=7, height=6, res=600)
  geo_time_plot
dev.off()

############################ quaSSE estimates ###################################

### loading estimates from the best-fit quasse model
lm_lin_params = read.table("6_quasse/lm_lin_params.csv", sep=",", h = T)
# setting the common parameters
common_params = lm_lin_params

### common values
x1 = rep(0, nrow(common_params))
xend = rep(max(sqrt(hvolumes) ), nrow(common_params))

###  mean line 
intercepts = common_params$l.c 
angulars = common_params$l.m
lines = data.frame(x1, xend, y1= intercepts, yend= xend*angulars + intercepts)
# mean line
mean_line = apply(lines, MARGIN = 2, FUN= mean)
# se line
ci_line = apply(lines, MARGIN = 2, FUN= sd) /  sqrt(nrow(lines)) * 1.96
ci_line_up = mean_line + ci_line
ci_line_dw =  mean_line - ci_line
# dataframe with all lines
lines_df = data.frame(rbind(mean_line, ci_line_up, ci_line_dw))
# mean function to get lambda
get_lambda = function(x){ (mean(angulars)*x) + mean(intercepts)}

### hypervolume and lambda estimate per group
mean_hv = aggregate(hvolumes, by=list(geo_states), mean)
ci_hv = aggregate(hvolumes, by=list(geo_states), function(x){ 1.96*(sd(x)/ sqrt(length(x))) })
## getting lambda estimates
mean_lambda = get_lambda(mean_hv$x)
ci_lambda = get_lambda(ci_hv$x)
lambda_estimates = data.frame(mean_lambda, ci_lambda)

### plotting
# axis title size
axis_title_size = 10
# axis names
x_axis_name = "hypervolume size"
y_axis_name = "lambda"

qua_plot = ggplot() +
  geom_segment(data = lines_df[1,], color="black", size=0.75, alpha=1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = lines_df[2,], color="black", size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_segment(data = lines_df[3,], color="black", size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  xlab(x_axis_name) + ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text=element_text(size=6),legend.position = "none")

tiff("7_graphs/quasse_estimates.tiff", units="cm", width=7, height=6, res=600)
  qua_plot
dev.off()

############################## SUPPLEMENTARY MATERIAL ######################

### loading fit values and model parameters
best_fit_models = read.table("geosse_time/geosse_time_best_fit_models.csv", sep=",", h = T)
sd_s_lin_params = read.table("geosse_time/sd_s_lin_params.csv", sep=",",  h = T)

### most common model
most_repeated = max(table(best_fit_models$model_name))
frequent_model = names(which(table(best_fit_models$model_name) == most_repeated) ) 
frequent_model_n = which(best_fit_models$model_name == frequent_model)
frequent_model_params = sd_s_lin_params[frequent_model_n,]

### common values
x1 = rep(0, nrow(frequent_model_params))
xend = rep(5, nrow(frequent_model_params))
y_lim=c(0,1.5)

# axis names
x_axis_name = "time (million years)"
y_axis_name = "lambda"

### sA line estimates = outside
intercepts_a = frequent_model_params$sA.c 
angulars_a = frequent_model_params$sA.m
lines_a = data.frame(x1, xend, y1= intercepts_a, yend= xend*angulars_a + intercepts_a)

sa = ggplot() +
  geom_segment(data = lines_a, color="#D81B60", size=1.2, alpha=0.1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  ylim(y_lim)+
  xlab(x_axis_name) + ylab("")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=9,face="bold"),axis.text=element_text(size=6),legend.position = "none")

### sB line estimates = AF
intercepts_b = frequent_model_params$sB.c 
angulars_b = frequent_model_params$sB.m
lines_b = data.frame(x1, xend, y1= intercepts_b, yend = xend*angulars_b + intercepts_b)

sb = ggplot() +
  geom_segment(data = lines_b, color="#1E88E5", size=1.2, alpha=0.1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  ylim(y_lim)+
  xlab(x_axis_name)+ ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=9,face="bold"),axis.text=element_text(size=6),legend.position = "none")

### sAB line estimates
intercepts_ab = common_params$sAB.c 
angulars_ab = common_params$sAB.m
lines_ab = data.frame(x1, xend, y1= intercepts_ab, yend = xend*angulars_ab + intercepts_ab)

sab = ggplot() +
  geom_segment(data = lines_ab, color="#FFC107", size=1.2, alpha=0.1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  ylim(y_lim)+
  xlab(x_axis_name)+ ylab("")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=9,face="bold"),axis.text=element_text(size=6),legend.position = "none")

tiff("geosse_time/geosse_time_lambda.tiff", units="in", width=4.5, height=1.75, res=600)
ggarrange(sb,sab, sa, nrow=1,ncol=3)
dev.off()

