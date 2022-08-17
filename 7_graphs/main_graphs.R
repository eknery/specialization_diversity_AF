setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

library (phytools)
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

################################## Phylo tree and LTT plot ####################################

### stochastic mapping
# simulate character evolution
maps = make.simmap(mcc, geo_states , model="ARD", pi='estimated', nsim=100)

# choos one map
n = 10

# count lineage through time
ltt_obj = ltt(maps[[n]], plot=FALSE)

# plot both!
layout(matrix(c(1,2),2,1),heights=c(0.5,0.4))
  plot(maps[[n]], mycols, ftype="off", mar=c(0,4.1,1.1,1.1))
  grid( ny= NA)
  par(mar=c(5.1,4.1,0,1.1))
  plot(ltt_obj,colors=mycols,bty="n",las=1,lwd=3,show.tree=F, legend=T, ylim=c(1,Ntip(mcc)), xlab="time (above the root)")
  legend("topleft",names(mycols),pch=22,pt.bg=mycols, bty="n",cex=0.8,pt.cex=1.2)
  grid( ny= NA)
simmap_ltt_plot = recordPlot()

tiff("7_graphs/simmap_ltt.tiff", units="in", width=4, height=6, res=600)
  simmap_ltt_plot
dev.off()

####################### Describing altitude, time, and niche breadth #############

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
  xlab("geographic distribution")+ ylab("median altitude (m)")+
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

tiff("7_graphs/species_data.tiff", units="in", width=3, height=5, res=600)
  ggarrange(alt_plot,hv_plot, nrow=2,ncol=1)
dev.off()

################################## RO and NO intercepts #########################

### loading geographic and niche sister data
sister_no_metrics = read.table("2_sister_hypervolume/sister_no_metrics.csv", sep=',', h=T)
sister_ro_metrics = read.table("3_sister_geography/sister_ro_metrics.csv", sep=',', h=T)

axis_title_size = 10

### sister ro metrics
ro_plot = ggplot(data= sister_ro_metrics, aes(x=state, y=intercept_ro, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.50) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("RO intercept")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")

### sister no metrics
no_plot = ggplot(data= sister_no_metrics, aes(x=state, y=intercept_no, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.50) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("NO intercept")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")

tiff("7_graphs/intercepts_data.tiff", units="in", width=6, height=2.5, res=600)
  ggarrange(ro_plot,no_plot, nrow=1,ncol=2)
dev.off()

####################### angular coefficients #################################

###### best-fit geosse-time estimates 

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
xend = rep(3.5, nrow(frequent_model_params))
y_lim=c(0,1)

# axis names
x_axis_name = "time (million years ago)"
y_axis_name = "lambda"

### only AF parameters

# sB line estimates = AF
intercepts_b = frequent_model_params$sB.c 
angulars_b = frequent_model_params$sB.m
lines_b = data.frame(x1, xend, y1= intercepts_b, yend = xend*angulars_b + intercepts_b)

# mean line
mean_line_b = apply(lines_b, MARGIN = 2, FUN= mean)

# se line
se_line_b = apply(lines_b, MARGIN = 2, FUN= sd) /  sqrt(nrow(lines_b))
se_line_b_up = mean_line_b + se_line_b
se_line_b_dw =  mean_line_b - se_line_b

# dataframe with all lines
lines_b_df = data.frame(rbind(mean_line_b, se_line_b_up, se_line_b_dw))

### plotting
axis_title_size = 10

geo_time_plot = ggplot() +
  geom_segment(data = lines_b_df[1,], color="black", size=0.75, alpha=1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = lines_b_df[2,], color="black", size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_segment(data = lines_b_df[3,], color="black", size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_vline(aes(xintercept = c(0.05,2.58)),linetype="dotted",colour="gray",size=0.75)+
  ylim(y_lim)+
  xlab(x_axis_name)+ ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size ,face="bold"),axis.text=element_text(size=6),legend.position = "none")

###### best-fit quasse estimates 

### loading estimates from the best-fit quasse model
lm_lin_params = read.table("quasse/lm_lin_params.csv", sep=",", h = T)

# setting the common parameters
common_params = lm_lin_params

### summary hypervolume per geographic state
mean_hv = aggregate(hvolumes, by=list(geo_states), mean)
sd_hv = aggregate(hvolumes, by=list(geo_states), sd)
n = aggregate(hvolumes, by=list(geo_states),length)
se_hv = sd_hv[,2]/sqrt(n[,2])
summary_hv = data.frame(mean_hv, sd_hv[,2], se_hv)
colnames(summary_hv)=c("state","mean","sd","se")

### common values
x1 = rep(0, nrow(common_params))
xend = rep(max(hvolumes), nrow(common_params))
y_lim=c(0,1)

### line estimates
intercepts = common_params$l.c 
angulars = common_params$l.m
lines = data.frame(x1, xend, y1= intercepts, yend= xend*angulars + intercepts)

# mean line
mean_line = apply(lines, MARGIN = 2, FUN= mean)

# se line
se_line = apply(lines, MARGIN = 2, FUN= sd) /  sqrt(nrow(lines))
se_line_up = mean_line + se_line
se_line_dw =  mean_line - se_line

# dataframe with all lines
lines_df = data.frame(rbind(mean_line, se_line_up, se_line_dw))

# axis names
x_axis_name = "hypervolume size"
y_axis_name = "lambda"

### plotting
axis_title_size = 10
qua_plot = ggplot() +
  geom_segment(data = lines_df[1,], color="black", size=0.75, alpha=1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = lines_df[2,], color="black", size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_segment(data = lines_df[3,], color="black", size=0.5, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  #geom_vline(data= summary_hv, aes(xintercept = mean), linetype="dotted", color= "#1E88E5", size=0.5)+
  geom_vline(data= summary_hv, aes(xintercept = mean-se),linetype="dotted", colour= c(mycols), size=0.5)+
  geom_vline(data= summary_hv, aes(xintercept = mean+se),linetype="dotted", colour= c(mycols), size=0.5)+
  ylim(y_lim)+
  xlab(x_axis_name) + ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=axis_title_size,face="bold"),axis.text=element_text(size=6),legend.position = "none")

tiff("7_graphs/angular_estimates.tiff", units="in", width=6, height=2.5, res=600)
  ggarrange(qua_plot,geo_time_plot, nrow=1,ncol=2)
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

