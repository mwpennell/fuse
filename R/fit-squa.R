## squamate analysis
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
## construct model
lik <- make.mkn(phy, states, k=4)
#con <- constrain(lik, q14~0, q24~0, q32~0, q42~0)

p <- starting.point.musse(phy, 4)
#p <- p[argnames(con)]

## mle 
#mle <- find.mle(con, p)

prior <- make.prior.exponential(10)

#tmp <- mcmc(con, x.init = p, w=1, prior=prior,
#            nsteps=100, print.every=0)

#w <- diff(sapply(tmp[2:9], range))

#samp <- mcmc(con, x.init = p, w=w, prior=prior,
#             nsteps=50000, print.every=1000)

#saveRDS(samp, "output/results/fish.rds")


## constrain even more
## assume transition from XY to ZW are equally likely regardless of whether one
## is fused or not

feq <- constrain(lik, q14~0, q24~0, q32~0, q42~0,
                 q43~qzw, q23~qzw, q13~qzw,
                 q41~qxy, q31~qxy, q21~qxy,
                 extra=c("qzw", "qxy"))
p.feq <- rep(p["q12"], 4)
names(p.feq) <- argnames(feq)

tmp <- mcmc(feq, x.init = p.feq, w=1, prior=prior,
            nsteps=100, print.every = 0)

w <- diff(sapply(tmp[2:5], range))

samp <- mcmc(feq, x.init = p.feq, w=1, prior=prior,
             nsteps=50000, print.every = 1000)

saveRDS(samp, "output/results/squa-feq.rds")


