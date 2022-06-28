setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

library(ape)
library(diversitree)
library(FAmle)

### loading phylogenetic tree
mcc = read.tree("0_data/mcc_phylo.nwk")
phylo_trees = read.tree("0_data/100_rand_phylos.nwk")

### loading trait dataset
spp_hvolumes = read.table("1_hypervolume_inference/spp_hvolumes.csv", h=T, sep=",")

### setting discrete states
summary(spp_hvolumes$hvolume)
ths_value = 2.69 
trait_states = spp_hvolumes$hvolume
trait_states[spp_hvolumes$hvolume <= ths_value] = 1 #specialist
trait_states[spp_hvolumes$hvolume > ths_value] = 0 #generalist
names(trait_states) = spp_hvolumes$species

table(trait_states)

### bisse
start_params = starting.point.bisse(mcc)

bisse_full = make.bisse(mcc,  trait_states)
bisse_m_const  =  constrain(bisse_full,  mu0  ~  mu1)
bisse_lm_const = constrain(bisse_full,  lambda0  ~  lambda1, mu0  ~  mu1)

mle_full = find.mle(bisse_full, start_params)
mle_m_const = find.mle(bisse_m_const, start_params)
mle_lm_const = find.mle(bisse_lm_const, start_params)


anova(mle_lm_const, mle_m_const, mle_full)



