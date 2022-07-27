setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

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

### loading geographic and niche sister data
sister_no_metrics = read.table("2_sister_hypervolume/sister_no_metrics.csv", sep=',', h=T)

### loading geographic and niche sister data
sister_ro_metrics = read.table("3_sister_geography/sister_ro_metrics.csv", sep=',', h=T)

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

# my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

### sister ro metrics
tiff("7_graphs/RO_intercept_geographic_distribution.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= sister_ro_metrics, aes(x=state, y=intercept_ro, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.50) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("RO intercept")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()

### sister no metrics
tiff("7_graphs/NO_intercept_geographic_distribution.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data= sister_no_metrics, aes(x=state, y=intercept_no, fill= state)) +
  geom_point(aes(color=state),position = position_jitter(width = 0.07), size = 1.5, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.50)+
  geom_flat_violin(position = position_nudge(x = 0.12, y = 0), alpha = 0.50) +
  scale_fill_manual(values=mycols)+
  scale_colour_manual(values=mycols)+
  xlab("geographic distribution")+ ylab("NO intercept")+
  scale_x_discrete(labels=c("AF" = "AF-endemic", "AFother" = "AF and other\ndomains", "other" = "outside AF"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=8),legend.position = "none")
dev.off()

############################## best-fit geosse-time estimates #################################

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

##### only AF-lambda

# mean line
mean_intercept_b = mean(intercepts_b)
mean_angular_b = mean(angulars_b)
mean_line_b = data.frame(x1, xend, y1= mean_intercept_b, yend= xend*mean_angular_b + mean_intercept_b)

# common values
x1 = rep(0, nrow(frequent_model_params))
xend = rep(5, nrow(frequent_model_params))
y_lim=c(0,1.5)

# axis names
x_axis_name = "time (million years)"
y_axis_name = "lambda"

tiff("7_graphs/geosse_time_lambda_AF.tiff", units="in", width=3.5, height=3, res=600)
ggplot() +
  geom_segment(data = lines_b, color="gray", size=0.75, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = mean_line_b, color="black", size=1, alpha=1,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_vline(aes(xintercept = 2.58),linetype="dotted",colour="black",size=0.75)+
  ylim(y_lim)+
  xlab(x_axis_name)+ ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=9,face="bold"),axis.text=element_text(size=6),legend.position = "none")
dev.off()


############################## best-fit quasse estimates #################################

### loading estimates from the best-fit quasse model
lm_lin_params = read.table("quasse/lm_lin_params.csv", sep=",", h = T)

# setting the common parameters
common_params = lm_lin_params

### loading species' hypervolumes
spp_hvolumes = read.table("1_hypervolume_inference/spp_hvolumes.csv", h=T, sep=",")

# named vector of hypervolumes
hvolumes = spp_hvolumes$hvolume
names(hvolumes) = spp_hvolumes$species

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

# define states threshold
high_ths = 0.9
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "AF"
geo_states[af_percentage <= low_ths] = "other"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "AFother"
names(geo_states) = spp_count_domain$species

# my colors
mycols = c( "#1E88E5", "#FFC107", "#D81B60")

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
y_lim=c(0,1.5)

### line estimates
intercepts = common_params$l.c 
angulars = common_params$l.m
lines = data.frame(x1, xend, y1= intercepts, yend= xend*angulars + intercepts)

mean_intercept = mean(intercepts)
mean_angular = mean(angulars)
mean_line = data.frame(x1, xend, y1= mean_intercept, yend= xend*mean_angular + mean_intercept)

# axis names
x_axis_name = "hypervolume size"
y_axis_name = "lambda"

tiff("7_graphs/quasse_lambda.tiff", units="in", width=3.5, height=3, res=600)
ggplot() +
  geom_segment(data = lines, color="gray", size=0.75, alpha=0.5,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  geom_segment(data = mean_line, color="black", size=1, alpha=1,  aes(x = x1, y = y1, xend = xend, yend = yend))+
  geom_vline(data= summary_hv[1,], aes(xintercept = mean),linetype="solid", color= "#1E88E5", size=0.75)+
  geom_vline(data= summary_hv[1,], aes(xintercept = mean-sd),linetype="dotted", color= "#1E88E5", size=0.5)+
  geom_vline(data= summary_hv[1,], aes(xintercept = mean+sd),linetype="dotted", color= "#1E88E5", size=0.5)+
  ylim(y_lim)+
  xlab(x_axis_name) + ylab(y_axis_name)+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=9,face="bold"),axis.text=element_text(size=6),legend.position = "none")
dev.off()
