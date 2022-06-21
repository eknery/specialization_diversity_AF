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


sister_hv_comparison = read.table("2_sister_hypervolume/sister_hv_comparison.csv", sep=",", h = T)
sister_hv_distance = read.table("2_sister_hypervolume/sister_hv_distance.csv", sep=',', h = T)

sister_area_comparison = read.table("3_sister_geography/sister_area_comparison.csv", sep=",", h = T)
sister_geo_disatnce = read.table("3_sister_geography/sister_geo_distance.csv", sep=",", h = T)



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

ggplot(data= sister_geo_distance, aes(x=-divergence_time, y=distance_to_sister, fill=state)) +
  geom_smooth(method= lm, formula= y~x, aes(color=state), fill="white")+
  geom_point(aes(color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))


### sympatry ~ hvolume
df = data.frame(sister_geo_distance, sister_hv_comparison$sp_hvolume)
colnames(df)[5] = "sp_hvolume"

ggplot(data= df, aes(x=sp_hvolume, y=distance_to_sister, fill=state)) +
  geom_smooth(method= lm, formula= y~x, aes(color=state, fill=state) )+
  geom_point(aes(color=state), position = position_jitter(width = 0.07), size = 2) +
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10))
