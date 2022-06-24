setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

library(ape)
library(diversitree)
library(phangorn)
library(FAmle)

source("traitDependent_functions.R")

### loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")

### loading trait dataset
spp_geographic_distribution = read.table("0_data/spp_geographic_distribution.csv", h=T, sep=",")

### setting geographic states
geo_states = spp_geographic_distribution$state
geo_states[geo_states == "AFother"] = 0
geo_states[geo_states == "other"] = 0
geo_states[geo_states == "AF"] = 1
geo_states = as.numeric(geo_states)
names(geo_states) = spp_geographic_distribution$species

res1 <- FISSE.binary(phy=mcc, states=geo_states, qratetype = "parsimony" )
pval_2tailed   <- min(res1$pval, 1-res1$pval)*2

### setting discrete categories
ths_value = 1.5 #median(spp_hvolumes$hvolume)
hv_states = spp_hvolumes$hvolume
hv_states[hv_states < ths_value] = "1" #specialist
hv_states[hv_states != "1"] = "0" #generalist
hv_states = as.numeric(hv_states)
names(hv_states) = spp_hvolumes$species

res2 <- FISSE.binary(phy=mcc, states=hv_states, qratetype = "parsimony" )
pval_2tailed   <- min(res2$pval, 1-res2$pval)*2
