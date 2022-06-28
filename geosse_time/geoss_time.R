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

### time functions
t_functions=rep(c("linear.t", "constant.t"),c(3, 4))

### set models
geosse_t_full = make.geosse.t(mcc, states=geo_states, functions=t_functions, sampling.f=0.9)
geosse_t_x_const = constrain(geosse_t_full, xA ~ xB)
geosse_t_sx_const = constrain(geosse_t_full, sA.m ~ sB.m, xA ~ xB)

### starting values = sA.c, sA.m, sB.c, sB.m, sAB.c, sAB.m, xA, xB, dA, dB
p = starting.point.geosse(mcc)
start_geosse_t = c(p[1], 1/2*p[1], p[2], 2*p[2], p[3], 1.5*p[3], p[4], p[5], p[6], p[7])
names(start_geosse_t) = argnames(geosse_t_full)

### find mle
mle_geosse_t_full = find.mle(geosse_t_full, start_geosse_t)
mle_geosse_t_x_const = find.mle(geosse_t_x_const, start_geosse_t)
mle_geosse_t_sx_const = find.mle(geosse_t_sx_const, start_geosse_t)

anova(mle_geosse_t_sx_const, mle_geosse_t_x_const, mle_geosse_t_full)
mle_geosse_t_sx_const$par
mle_geosse_t_x_const$par

