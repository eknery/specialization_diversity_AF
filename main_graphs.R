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

### sympatry ~ hvolume
spp_geo_hv = data.frame(sister_geo_distance, sister_hv_comparison$sp_hvolume)
colnames(spp_geo_hv )[5] = "sp_hvolume"

ggplot(data= spp_geo_hv , aes(x=sp_hvolume, y=log(distance_to_sister) ) ) +
  geom_smooth(method= lm, formula= y~x )+
  geom_point(aes(color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))
