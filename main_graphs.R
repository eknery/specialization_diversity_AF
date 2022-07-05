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

### loading niche and geographic sister data
sister_hv_comparison = read.table("2_sister_hypervolume/sister_hv_comparison.csv", sep=",", h = T)
sister_geo_distance = read.table("3_sister_geography/sister_geo_distance.csv", sep=",", h = T)

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

sister_hv_comparison$state = geo_states
sister_geo_distance$state = geo_states

### hvolume overlap
ggplot(data= sister_hv_comparison, aes(x=state, y=intersection/union, fill= state)) +
  geom_boxplot(alpha=0.5)+
  geom_point(aes(y=intersection/union, color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))

ggplot(data= sister_hv_comparison, aes(x=divergence_time, y=intersection/union, fill=state)) +
  geom_smooth(method= lm, formula= y~x, aes(color=state), fill="white")+
  geom_point(aes(color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))

### geo distance !!!!!

ggplot(data= sister_geo_distance, aes(x=state, y=distance_to_sister, fill= state)) +
  geom_boxplot(alpha=0.5)+
  geom_point(aes(y=distance_to_sister, color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))

ggplot(data= sister_geo_distance, aes(x=divergence_time, y=distance_to_sister, fill=state)) +
  geom_smooth(method= lm, formula= y~x, aes(color=state), fill="white")+
  geom_point(aes(color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))


### sympatry ~ hvolume
spp_geo_hv = data.frame(sister_geo_distance, sister_hv_comparison$sp_hvolume)
colnames(spp_geo_hv )[5] = "sp_hvolume"

ggplot(data= spp_geo_hv , aes(x=sp_hvolume, y=log(distance_to_sister) ) ) +
  geom_smooth(method= lm, formula= y~x )+
  geom_point(aes(color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))
