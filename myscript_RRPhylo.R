install.packages("RRphylo")

### packages
library(ape)
library(phytools)
library(geiger)
library(OUwie)
library(RRphylo)

### loading scale  niche
hv_scale_positions = read.table("hv_scale_positions.csv", header =T, sep=",",  na.strings = "NA", fill=T)
scale_traits = hv_scale_positions[,2:5]
rownames(scale_traits) = hv_scale_positions[,1]

### loading hv volume
hv_spp_volumes = read.table("hv_spp_volumes.csv", header =T, sep=",",  na.strings = "NA", fill=T)
volumes = hv_spp_volumes$volumes
names(volumes) = hv_spp_volumes$species

### loading trees
mcc=read.tree("mcc_phylo.nwk")
trees=read.tree("100_rand_phylos.nwk")

### phylogenetic ridge regression
RR = RRphylo(tree=mcc,y=scale_traits)

multi_rates = RR$multiple.rates
rates = multi_rates[,3]
pres_rates = rates[names(rates) %in% mcc$tip.label]

distribution = names(pres_rates)
distribution[distribution %in% AF] ="AF-endemic"
distribution[distribution != "AF-endemic"] = "widespread"
names(distribution)= names(pres_rates)
rates_df= data.frame(distribution,pres_rates)

ggplot(data= rates_df, aes(x=distribution, y=pres_rates, fill=distribution)) +
  geom_point(aes(y=pres_rates, color=distribution), position = position_jitter(width = 0.07), size = 1, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.3)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")

##########################################################


phenogram(mcc, volumes ,spread.labels=TRUE,spread.cost=c(1,0))

########################################################
distribution = hv_scale_positions$species
distribution[distribution %in% AF] ="AF-endemic"
distribution[distribution != "AF-endemic"] = "widespread"

hv_scale_positions = data.frame(distribution, hv_scale_positions)

ggplot(data= hv_scale_positions, aes(x=temperature_range, y=precipitation_seasonality, fill=distribution)) +
  geom_point(aes(color=distribution), size = 2, alpha = 0.75) +
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")


