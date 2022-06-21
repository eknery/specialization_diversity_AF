setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

install.packages("FAmle")

library(diversitree)
library(FAmle)

### loading phylogenetic tree
mcc = read.tree("mcc_phylo.nwk")

### loading trait dataset
hv_spp_volumes = read.table("hv_spp_volumes.csv", h=T, sep=",")

### setting trait vector
hvolume = hv_spp_volumes$volumes
names(hvolume) = hv_spp_volumes$species

### setting distribution vector
#defining endemics
AF=c("atlantica","baumgratziana","brunnea","budlejoides","capixaba", 
     "castaneiflora","cinerascens","discolor","dura","fasciculata",
     "formosa","hyemalis","kollmannii","kriegeriana","lymanii","mellina",
     "octopetala","penduliflora","petroniana","polyandra","racemifera",
     "robusta","ruschiana","setosociliata","shepherdii","valtheri","willdenowii")

distribution = hv_spp_volumes$species
distribution[hv_spp_volumes$species %in% AF] ="AF-endemic"
distribution[!hv_spp_volumes$species %in% AF] ="widespread"

### testing trait normality
hist(sqrt(hvolume))
shapiro.test(sqrt(hvolume))
sq_hvolume = sqrt(hvolume)

### sampling error
sq_hvolume_se = sd(sq_hvolume)/length(sq_hvolume)

### starting parameter values
p = starting.point.quasse(mcc, sq_hvolume)
p

###
exp_effect = max(sq_hvolume) - min(sq_hvolume) / min(sq_hvolume)

### setting linear function
xr = c(0, (max(sq_hvolume) + exp_effect) )
linear.x = make.linear.x(x0=xr[1], x1=xr[2])

### setting quasse functions
quasse_const =  make.quasse(mcc, sq_hvolume, sq_hvolume_se, lambda=constant.x, mu=constant.x, sampling.f=0.8)
quasse_linear = make.quasse(mcc, sq_hvolume, sq_hvolume_se, lambda=linear.x, mu=constant.x, sampling.f=0.8)
quasse_sigm = make.quasse(mcc, sq_hvolume, sq_hvolume_se, lambda=sigmoid.x, mu=constant.x, sampling.f=0.8)

############################# maximum likelihood estimates ##################

### constant function
argnames(quasse_const)
init_const = c(p[1],p[2],0.01, p[3])
names(init_const) = argnames(quasse_const)
lower_const = c(0,0,0,0)
names(lower_const) = argnames(quasse_const)
# finding constant mle
const_mle = find.mle(quasse_const, x.init=init_const, lower=lower_const)
const_mle$par

### linear function
argnames(quasse_linear)
init_linear = c(p[1], exp_effect*p[1], p[2], 0.01, p[3])
names(init_linear)= argnames(quasse_linear)
lower_linear = c(0,0,0,0,0)
names(lower_linear) = argnames(quasse_linear)
# finding linear mle
linear_mle = find.mle(quasse_linear, x.init=init_linear, lower=lower_linear)
linear_mle$par

### sigmoid function
argnames(quasse_sigm)
init_sigm = c(p[1], exp_effect*p[1], mean(sq_hvolume),0.5, p[2],0.01, p[3])
names(init_sigm)= argnames(quasse_sigm)
lower_sigm = c(0,0,0,0,0,0,0)
names(lower_sigm) = argnames(quasse_sigm)
# finding linear mle
sigm_mle = find.mle(quasse_sigm, x.init=init_sigm, lower=lower_sigm)
sigm_mle$par

anova(const_mle, linear_mle, sigm_mle)

################################ plotting ####################################
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

a = linear_mle$par[2] - linear_mle$par[1]/xr[2]- xr[1]
b = linear_mle$par[1] 

lambda_f = function(x){a*x + b}
df = data.frame(distribution, sq_hvolume, lambda=lambda_f(sq_hvolume))

ggplot(data= df, aes(y=lambda, x=sq_hvolume, fill=distribution)) +
  geom_point(aes(color=distribution), position = position_jitter(width =linear_mle$par[5] ), size = 3, alpha = 0.25) +
  geom_abline(intercept = b, slope = a, color="gray", linetype="dashed", size=1.5)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  labs(x="sqrt(hypervolume)", y="lambda")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=10),legend.position = "none")

ggplot(data= df, aes(x=distribution,y=lambda, fill=distribution)) +
  geom_point(aes(y = lambda, color =distribution), position = position_jitter(width = 0.07), size = 2.5, alpha = 0.5) +
  geom_flat_violin(col="black", position = position_nudge(x = 0.12, y = 0), alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.25)+
  scale_fill_manual(values=c("green3","magenta"))+
  scale_colour_manual(values=c("green3","magenta"))+
  ylim(c(0.15,0.4))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12),legend.position = "none")
