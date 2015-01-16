## fish analysis
library(diversitree)

td <- readRDS("output/data/mamm.rds")

phy <- td$phy

## remove polytomies
phy <- multi2di(phy)

dat <- td$data

lik <- make.mkn(phy, dat, k=3)

p <- starting.point.musse(phy, 3)

lik <- constrain(lik, q23~0, q32~0)

p <- p[argnames(lik)]

prior <- make.prior.exponential(10)

tmp <- mcmc(lik, x.init = p, w=1, prior=prior,
            nsteps=100, print.every=0)

w <- diff(sapply(tmp[2:5], range))

samp <- mcmc(lik, x.init = p, w=w, prior=prior,
             nsteps=50000, print.every=1000)

saveRDS(samp, "output/results/mamm.rds")
