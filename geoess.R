setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

library(ape)
library(diversitree)
library(FAmle)

### loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

### loading trait dataset
spp_geographic_distribution = read.table("0_data/spp_geographic_distribution.csv", h=T, sep=",")

### setting geoess states
geo_states = spp_geographic_distribution$state
geo_states[geo_states == "AFother"] = 0
geo_states[geo_states == "other"] = 1
geo_states[geo_states == "AF"] = 2
geo_states = as.numeric(geo_states)
names(geo_states) = spp_geographic_distribution$species

### 
start_geosse = starting.point.geosse(mcc)

geosse_full = make.geosse(mcc, states=geo_states, sampling.f=0.9)
geosse_x_const = constrain(geosse_full, xA ~ xB, sampling.f=0.9)
geosse_sx_const = constrain(geosse_full, sA ~ sB, xA ~ xB, sampling.f=0.9)


mle_geosse_full = find.mle(geosse_full, start_geosse)
mle_geosse_x_const = find.mle(geosse_x_const, start_geosse)
mle_geosse_sx_const = find.mle(geosse_sx_const, start_geosse)

anova(mle_geosse_sx_const, mle_geosse_x_const, mle_geosse_full)
