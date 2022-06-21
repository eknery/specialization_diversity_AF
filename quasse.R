library(diversitree)

lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5) 
mu<- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)


function(x) linear.x(x, 0.1, 0.2, 0, 2.5)

phy <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0, single.lineage=FALSE)
str(phy)
phy$orig

states <- phy$tip.state 
states.sd <- 1/200

lik <- make.quasse(phy, states, states.sd, sigmoid.x, constant.x)

p<- starting.point.quasse(phy, states) 
p

lik.nodrift <- constrain(lik, drift ~ 0) 
argnames(lik)
argnames(lik.nodrift)

p.start <- c(p[1], p[1], mean(states), 1, p[2:3]) 
names(p.start) <- argnames(lik.nodrift) 
p.start

lower <- c(0, 0, min(states),-Inf, 0, 0)

fit <- find.mle(lik.nodrift, p.start, control=list(parscale=.1), lower=lower, verbose=0)


