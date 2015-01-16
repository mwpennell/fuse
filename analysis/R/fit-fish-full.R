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

## 1 is XY unfused
## 2 is XY fused
## 3 is ZW unfused
## 4 is ZW fused

## construct model
lik <- make.mkn(phy, states, k=4)

p <- starting.point.musse(phy, 4)

## 1. Disallow ZW -> XXY and XY -> ZZW and XXY <-> ZZW
lik.con1 <- constrain(lik, q14~0, q24~0, q32~0, q42 ~ 0)
p.con1 <- p[argnames(lik.con1)]
#m1 <- find.mle(lik.con1, p.con1)

## 2. Constrain ZZW -> XY to equal ZW -> XY; same with XXY -> ZW
lik.con2 <- constrain(lik.con1, q41~q31, q23~q13)
p.con2 <- p[argnames(lik.con2)]
#m2 <- find.mle(lik.con2, p.con2)

lik.con2a <- constrain(lik.con1, q41~0, q23~0)
p.con2a <- p[argnames(lik.con2a)]
#m2a <- find.mle(lik.con2a, p.con2a)

## 3. Constrain XXY -> XY to be equal to ZZW -> ZW
lik.con3 <- constrain(lik.con2, q43~q21)
p.con3 <- p[argnames(lik.con3)]
#m3 <- find.mle(lik.con3, p.con3)

lik.con3a <- constrain(lik.con2, q43~0)
#p.con3a <- m1$par[argnames(lik.con3a)]
#m3a <- find.mle(lik.con3a, p.con3a)

lik.con3b <- constrain(lik.con2, q21~0)
#p.con3b <- m1$par[argnames(lik.con3b)]
#m3b <- find.mle(lik.con3b, p.con3b)

## 4. Are XY -> XXY different from ZW -> ZWW
lik.con4 <- constrain(lik.con3, q34~q12)
p.con4 <- p[argnames(lik.con4)]
#m4 <- find.mle(lik.con4, p.con4)

prior <- make.prior.exponential(20)

tmp <- mcmc(lik.con3, x.init = p.con3, w=1, prior=prior,
            nsteps=100, print.every=0)

w <- diff(sapply(tmp[2:6], range))

samp <- mcmc(lik.con3, x.init = p.con3, w=w, prior=prior,
             nsteps=50000, print.every=1000)

saveRDS(samp, "output/results/fish-final.rds")
