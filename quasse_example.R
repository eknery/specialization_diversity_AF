phy <- read.nexus("Vos-2006.nex") 
d<- read.csv("Redding-2010.csv") 
mass <- log(d$mass) 
names(mass) <- d$tip.label

mass.sd <- 1/50

p<- starting.point.quasse(phy, mass) 
p

xr<- range(mass) + c(-1,1) * 20 * p["diffusion"]

linear.x = make.linear.x(xr[1], xr[2])

make.primates <- function(lambda, mu){
  make.quasse(phy, mass, mass.sd, lambda, mu)
}

nodrift <- function(f){
  constrain(f, drift ~ 0)
}

f.c <- make.primates(constant.x, constant.x) 
f.l <- make.primates(linear.x, constant.x) 
f.s <- make.primates(sigmoid.x, constant.x)


control <- list(parscale=.1, reltol=0.001)
