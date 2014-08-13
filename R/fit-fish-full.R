## fish analysis
library(diversitree)

td <- readRDS("output/data/fish.rds")

phy <- td$phy
dat <- data.frame(td$data)

## reformat data for using regular MK models
states <- sapply(seq_len(nrow(dat)), function(x) paste(dat[x,1], dat[x,2], sep=""))
states <- factor(states, levels=c("00","01","10","11"), labels=c(1:4))
states <- as.numeric(states)
names(states) <- rownames(dat)

## construct model
lik <- make.mkn(phy, states, k=4)

p <- starting.point.musse(phy, 4)
p <- p[argnames(lik)]
prior <- make.prior.exponential(20)

tmp <- mcmc(lik, x.init = p, w=1, prior=prior,
            nsteps=100, print.every=0)

w <- diff(sapply(tmp[2:13], range))

samp <- mcmc(lik, x.init = p, w=w, prior=prior,
             nsteps=50000, print.every=1000)

saveRDS(samp, "output/results/fish-full.rds")
