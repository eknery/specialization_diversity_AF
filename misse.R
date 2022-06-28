# http://speciationextinction.info/articles/MiSSE-vignette.html

library(hisse)

turnover <- c(1,2)
eps <- c(0,0)
two_rate <- MiSSE(mcc, f=1, turnover=turnover, eps=eps, fixed.eps=0.9)
