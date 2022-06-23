setwd("C:/Users/eduar/Documents/GitHub/specialization_diversity_AF")

### packages
library(ape)
library(diversitree)

make.bd.t(tree, functions, sampling.f=NULL, unresolved=NULL,
          control=list(), truncate=FALSE, spline.data=NULL)


phy=read.tree("0_data/mcc_phylo.nwk")


pars <- c(0.1, 0.03)

lik.direct <- make.bd(phy)
lik.ode <- make.bd(phy, control=list(method="ode"))
lik.t <- make.bd.t(phy, c("constant.t", "constant.t"))

lik.direct(pars) - lik.ode(pars)


lik.t2 <- make.bd.t(phy, c("linear.t", "constant.t"))

pars2 <- c(pars[1], 0, pars[2])
lik.t2(pars2) - lik.t(pars)

fit <- find.mle(lik.direct, pars)
fit.t2 <- find.mle(lik.t2, pars2)

anova(fit, time.varying=fit.t2)
