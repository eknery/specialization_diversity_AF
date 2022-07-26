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

d_lin_params = read.table("geosse_time/d_lin_params.csv", sep=",",  h = T)
sd_s_lin_params = read.table("geosse_time/sd_s_lin_params.csv", sep=",",  h = T)
sxd_sx_lin_params = read.table("geosse_time/sxd_sx_lin_params.csv", sep=",", h = T)

### most common model
most_repeated = max(table(best_fit_models$model_name))
common_model = names(which(table(best_fit_models$model_name) == most_repeated) ) 
common_scenario = which(best_fit_models$model_name == common_model)
common_params = sd_s_lin_params[common_scenario,]

### common values
x1 = rep(0, nrow(common_params))
xend = rep(5, nrow(common_params))
y_lim=c(0,2)

# axis names
x_axis_name = "divergence time\n(million years)"
y_axis_name = "lambda"

### sA line estimates = outside
intercepts_a = common_params$sA.c 
angulars_a = common_params$sA.m
lines_a = data.frame(x1, xend, y1= intercepts_a, yend= xend*angulars_a + intercepts_a)

sa = ggplot() +
  geom_segment(data = lines_a, color="#D81B60", size=1.2, alpha=0.1,  aes(x = x1, y = y1, xend = xend, yend = yend), inherit.aes = FALSE)+
  ylim(y_lim)+
  xlab(x_axis_name) + ylab("")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=9,face="bold"),axis.text=element_text(size=6),legend.position = "none")

### sB line estimates = AF
intercepts_b = common_params$sB.c 
angulars_b = common_params$sB.m
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

############################## best-fit quasse estimates #################################

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


