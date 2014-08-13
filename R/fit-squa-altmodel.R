## squamate analysis (futher constraints)
library(diversitree)

td <- readRDS("output/data/squa.rds")

phy <- td$phy
dat <- data.frame(td$data)

## reformat data for using regular MK models
states <- sapply(seq_len(nrow(dat)), function(x) paste(dat[x,1], dat[x,2], sep=""))
states <- factor(states, levels=c("00","01","10","11"), labels=c(1:4))
states <- as.numeric(states)
names(states) <- rownames(dat)

## construct model
lik <- make.mkn(phy, states, k=4)

p <- starting.point.musse(phy,k=4)


prior <- make.prior.exponential(10)

feq <- constrain(lik, q14~0, q24~0, q32~0, q42~0,
                 q43~0, q23~0,
                 q41~0)
p.feq <- rep(p["q12"], 5)
names(p.feq) <- argnames(feq)

tmp <- mcmc(feq, x.init = p.feq, w=1, prior=prior,
            nsteps=100, print.every = 0)

w <- diff(sapply(tmp[2:6], range))

samp <- mcmc(feq, x.init = p.feq, w=1, prior=prior,
             nsteps=50000, print.every = 1000)

saveRDS(samp, "output/results/squa-5par.rds")
