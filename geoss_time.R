setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

library(ape)
library(diversitree)
library(FAmle)

### loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

### loading occurrence count per domain
spp_count_domain = read.table("0_data/spp_count_domain.csv", h=T, sep=",")

### define states threshold
high_ths = 0.75
low_ths = (1 - high_ths)
geo_states = af_percentage = spp_count_domain$AF/ apply(spp_count_domain[,-1], MARGIN = 1, FUN=sum)
geo_states[af_percentage >= high_ths] = "2"
geo_states[af_percentage <= low_ths] = "1"
geo_states[af_percentage > low_ths & af_percentage < high_ths] = "0"
geo_states = as.numeric(geo_states)
names(geo_states) = spp_count_domain$species

### 
start_geosse = starting.point.geosse(mcc)

geosse_full = make.geosse(mcc, states=geo_states, sampling.f=0.9)
geosse_x_const = constrain(geosse_full, xA ~ xB)
geosse_sx_const = constrain(geosse_full, sA ~ sB, xA ~ xB)


mle_geosse_full = find.mle(geosse_full, start_geosse)
mle_geosse_x_const = find.mle(geosse_x_const, start_geosse)
mle_geosse_sx_const = find.mle(geosse_sx_const, start_geosse)

anova(mle_geosse_sx_const, mle_geosse_x_const, mle_geosse_full)
mle_geosse_x_const$par
mle_geosse_full$par