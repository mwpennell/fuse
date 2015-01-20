## # Phylogenetic analyses of sex chromosome fusions

library(diversitree)
source("R/fusion-utility.R")

knitr::opts_chunk$set(tidy=FALSE)


## set colors
col.tree <- list(G=c("#86B8B1","#FA2A00"), K=c("white","#3D1C00"))
col.vlin <- "#3D1C00"
col.taxa <- c("#86B8B1", "#FA2A00")
names(col.taxa) <- c("fish", "squa")


## ## Part A: Counting fusions from raw data
## Read in all raw data
raw_all <- raw_karyotype_dat()

## Subset by taxonomic group
raw_fish <- subset(raw_all, Higher.taxonomic.group=="Fish")
raw_squa <- subset(raw_all, Order=="Squamata")
raw_chel <- subset(raw_all, Order=="Chelonia")
raw_mamm <- subset(raw_all, Higher.taxonomic.group=="Mammalia")
raw_aves <- subset(raw_all, Higher.taxonomic.group=="Aves")
raw_amph <- subset(raw_all, Higher.taxonomic.group=="Amphibia")


get_fusion_counts <- function(x){
   ct <- lapply(c("X2", "Y2", "Z2", "W2"), function(i) count_fusions(x,i))
   ct <- do.call(cbind, ct)
   xy_sys <- x[x$karyotype == "XY" | x$karyotype == "complex XY", ]
   xy_rec <- nrow(xy_sys)
   xy_uni <- length(unique(xy_sys$binom))
   zw_sys <- x[x$karyotype == "ZW" | x$karyotype == "complex ZW", ]
   zw_rec <- nrow(zw_sys)
   zw_uni <- length(unique(zw_sys$binom))
   tot <- nrow(x[x$karyotype != "",])
   uni <- length(unique(x[x$karyotype != "","binom"]))
   out <- cbind(ct, xy_rec, xy_uni, zw_rec, zw_uni, tot, uni)
   colnames(out) <- c("Y-A", "X-A", "W-A", "Z-A",
                      "XY_records", "XY_unique_spp",
                      "ZW_records", "ZW_unique_spp",
                      "total_records", "total_unique_spp")
   out
}

count_fusions <- function(x,i){
    id <- grep(i, x$notes)
    nm <- sapply(x$binom[id], function(i) {strsplit(i, split=" ")[[1]][1]})
    length(unique(nm))
}

raw_comb <- list(raw_fish, raw_squa, raw_chel, raw_mamm, raw_aves, raw_amph)

counts <- lapply(raw_comb, get_fusion_counts)
counts <- do.call(rbind, counts)
rownames(counts) <- c("Fish", "Squamates", "Turtles", "Mammals",
                      "Birds", "Amphibians")

write.csv(counts, "output/results/counts-raw.csv")


## ## Part B: XY vs ZW Systems

## define constraint functions

## 1. Disallow ZW -> XXY and XY -> ZZW and XXY <-> ZZW
constrain_doub <- function(lik)
    constrain(lik, q14~0, q24~0, q32~0, q42 ~ 0)

## 2. Constrain ZZW -> XY to equal ZW -> XY; same with XXY -> ZW
constrain_tnvr <- function(lik)
    constrain(lik, q41~q31, q23~q13)

## 3. Constrain XXY -> XY to be equal to ZZW -> ZW
constrain_back <- function(lik)
    constrain(lik, q43~q21)

## 4. Constrain XY fusions to occur at the same rate as ZW fusions
constrain_fuse <- function(lik)
    constrain(lik, q34~q12)

## find.mle with starting point (wrapper function)
fit_mle <- function(lik, sp)
    find.mle(lik, sp[argnames(lik)])

## mcmc with starting point (wrapper function)
fit_bay <- function(lik, sp, path){
    prior <- make.prior.exponential(20)
    tmp <- mcmc(lik, x.init=sp[argnames(lik)], w=1, prior=prior,
                nsteps=100, print.every=0)
    np <- length(argnames(lik))
    w <- diff(sapply(tmp[2:(np+1)], range))
    res <- mcmc(lik, x.init=sp[argnames(lik)], w=w, prior=prior,
                nsteps=50000, print.every=1000)
    saveRDS(res, path)
}


## Read in data
f_dat <- fish_karyotype_dat()
s_dat <- squa_karyotype_dat()

## Plot distribution on tree

## Note: This is slightly different from the figure as it appears in the main text
## Need to adjust this
fig_phylogeny <- function(dat){
    genus <- sub("_.+$", "", dat$phy$tip.label)
    trait.plot(ladderize(dat$phy), dat=as.data.frame(dat$data), cols=col.tree,
           lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")),
           class=genus, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)
}


fig_phylogeny(f_dat)
fig_phylogeny(s_dat)





## fish data
f_kary <- process_karyotype(f_dat)

## squamate data
s_kary <- process_karyotype(s_dat)



## ## Fish karyotype analyses (ML)

## define starting points
f_sp <- starting.point.musse(f_kary$phy, k=4)
## build likelihood function
f_lik <- make.mkn(f_kary$phy, f_kary$states, k=4)


## Add constraints in piecemeal
f_lik_doub <- constrain_doub(f_lik)
f_res_doub <- fit_mle(f_lik_doub, f_sp)

f_lik_tnvr <- constrain_tnvr(f_lik_doub)
f_res_tnvr <- fit_mle(f_lik_tnvr, f_sp)

anova(f_res_doub, f_res_tnvr)
## p=0.54
## Non-significant

f_lik_back <- constrain_back(f_lik_tnvr)
f_res_back <- fit_mle(f_lik_back, f_sp)

anova(f_res_tnvr, f_res_back)
## p=0.80
## Non-significant

f_lik_fuse <- constrain_fuse(f_lik_back)
f_res_fuse <- fit_mle(f_lik_fuse, f_sp)

anova(f_res_back, f_res_fuse)
## p=0.014
## Significant, fusion rates are different


## ## Fish karyotype analyses (Bayesian)

## fit f_lik_back
## fit_bay(f_lik_back, f_sp, "output/results/karyotype-fish.rds")


## Read in fusion chain
f_mcmc <- readRDS("output/results/karyotype-fish.rds")
## Remove burnin
f_mcmc <- f_mcmc[-seq_len(10000),]

## Differences between XY and ZW fusion rates
f_mcmc$diff <- f_mcmc$q12 - f_mcmc$q34

fig_fusion_rates <- function(x, taxa){
    profiles.plot(x["diff"], col.line=col.taxa[taxa], opacity=0.9,
                  xlab="Difference between XY and ZW fusion rates",
                  frame.plot=FALSE)
    axis(side=1)
    abline(v=0, lwd=2, lty=2, col=col.vlin)
}

fig_fusion_rates(f_mcmc, "fish")

## What proportion of posterior supports XY > ZW
length(which(f_mcmc$diff > 0))/nrow(f_mcmc)
## 98.4% of posterior probability
## 95% CI does not overlap with 0





## ## Squamate karyotype analyses (ML)

## define starting points
s_sp <- starting.point.musse(s_kary$phy,k=4)
## build likelihood function
s_lik <- make.mkn(s_kary$phy, s_kary$states, k=4)


## Add constraints in piecemeal
s_lik_doub <- constrain_doub(s_lik)
#s_res_doub <- fit_mle(s_lik_doub, s_sp)

s_lik_tnvr <- constrain_tnvr(s_lik_doub)
s_res_tnvr <- fit_mle(s_lik_tnvr, s_sp)

#anova(s_res_doub, s_res_tnvr)

s_lik_back <- constrain_back(s_lik_tnvr)
s_res_back <- fit_mle(s_lik_back, s_sp)

anova(s_res_tnvr, s_res_back)
## p=0.012
## This is significant due to the back transtion from ZZW -> ZW (q43) being v. high
s_res_tnvr$par

s_lik_fuse <- constrain_fuse(s_lik_back)
s_res_fuse <- fit_mle(s_lik_fuse, s_sp)

anova(s_res_back, s_res_fuse)
## p=0.003
## This is also highly signficant

## ## Squamate karyotype analyses (Bayesian)
## fit s_lik_back and s_lik_tnvr
# fit_bay(s_lik_back, s_sp, "output/results/karyotype-squa.rds")
# fit_bay(s_lik_tnvr, s_sp, "output/results/karyotype-squa-6par.rds")


## Two versions of analysis

## First, not estimating back transtion parameters separately

## Read in fusion chain
s_mcmc <- readRDS("output/results/karyotype-squa.rds")
## Remove burnin
s_mcmc <- s_mcmc[-seq_len(10000),]

## Differences between XY and ZW fusion rates
s_mcmc$diff <- s_mcmc$q12 - s_mcmc$q34

fig_fusion_rates(s_mcmc, "squa")

## What proportion of posterior supports XY > ZW
length(which(s_mcmc$diff > 0))/nrow(s_mcmc)
## 99.9% of posterior probability supports this conclusion



## Second, estimating back transtion parameters

## Read in fusion chain
s_mcmc_rev <- readRDS("output/results/karyotype-squa-6par.rds")
## Remove burnin
s_mcmc_rev <- s_mcmc_rev[-seq_len(10000),]

## Differences between XY and ZW fusion rates
s_mcmc_rev$diff <- s_mcmc_rev$q12 - s_mcmc_rev$q34

fig_fusion_rates(s_mcmc_rev, "squa")

## What proportion of posterior supports XY > ZW
length(which(s_mcmc_rev$diff > 0))/nrow(s_mcmc_rev)
## 92.0% of posterior probability supports XY > ZW
## This is not quite significant, owing to the high estimate for q43

## Here we also calculate the residency time
s_mcmc_rev$xy <- s_mcmc_rev$q12 / (s_mcmc_rev$q12 + s_mcmc_rev$q21)
s_mcmc_rev$zw <- s_mcmc_rev$q34 / (s_mcmc_rev$q34 + s_mcmc_rev$q43)
s_mcmc_rev$resid <- s_mcmc_rev$xy - s_mcmc_rev$zw

fig_fusion_resid <- function(x, taxa){
    profiles.plot(x["resid"], col.line=col.taxa[taxa], opacity=0.9,
                  xlab="Difference in fusion residency time (XY - ZW)",
                  frame.plot=FALSE)
    axis(side=1)
    abline(v=0, lwd=2, lty=2, col=col.vlin)
}

fig_fusion_resid(s_mcmc_rev, "squa")

## What proportion of posterior favors residency time of XY fusions
## being greater than residency time of ZW fusions?

length(which(s_mcmc_rev$resid > 0))/nrow(s_mcmc_rev)
## 99.8% of posterior probability. 95CI does not overlap








## ## PART C: Examining fusions on individual chromosomes

## read in individual chromosome data
f_chr <- fish_chromosome_dat()
s_chr <- squa_chromosome_dat()

## ## Fish chromosome analysis (ML)

## ### Definitions of states 
## 1: XY
## 2: ZW
## 3: YA (XXY)
## 4: XA (XYY)
## 5: ZA (ZWW)

## Note that ZZW is not present

## Baseline constraints (biologically enforced)

## 1. Disallow transtions from all fused to other fused states
## 2. Disallow transitions from XY <-> ZW fused / ZW <-> XY fused
constrain_chr_fish <- function(lik)
    constrain(lik, q34~0, q35~0, q43~0, q45~0, q53~0, q54~0,
              q15~0, q51~0, q23~0, q32~0, q24~0, q42~0)

## 3. Constrain rev transtions to all have equal rates
constrain_rev_fish <- function(lik)
    constrain(lik, q41~q31, q52~q31)

## 4a. Constrain YA fusions to be equal to XA
constrain_xyeq_fish <- function(lik)
    constrain(lik, q14~q13)

## 4b. Constrain YA fusions to be equal to ZA
constrain_yzeq_fish <- function(lik)
    constrain(lik, q25~q13)

## 4c. Constrain ZA fusions to be equal to XA
constrain_xzeq_fish <- function(lik)
    constrain(lik, q25~q14)


## define starting points
f_sp_chr <- starting.point.musse(f_chr$phy, k=5)

## build likelihood function
f_lik_chr <- make.mkn(f_chr$phy, f_chr$data, k=5)

## Baseline
f_lik_chr <- constrain_chr_fish(f_lik_chr)
f_res_chr <- fit_mle(f_lik_chr, f_sp_chr)

## Rev transitions
f_lik_rev <- constrain_rev_fish(f_lik_chr)
f_res_rev <- fit_mle(f_lik_rev, f_sp_chr)

anova(f_res_rev, f_res_chr)

## Keep this constraint for further analyses

## YA different from XA?
f_lik_xyeq <- constrain_xyeq_fish(f_lik_rev)
f_res_xyeq <- fit_mle(f_lik_xyeq, f_sp_chr)

anova(f_res_xyeq, f_res_rev)

## p=0.016*

## YA different from ZA?
f_lik_yzeq <- constrain_yzeq_fish(f_lik_rev)
f_res_yzeq <- fit_mle(f_lik_yzeq, f_sp_chr)

anova(f_res_yzeq, f_res_rev)

## p=0.035*

## XA different from ZA?
f_lik_xzeq <- constrain_xzeq_fish(f_lik_rev)
f_res_xzeq <- fit_mle(f_lik_xzeq, f_sp_chr)

anova(f_res_xzeq, f_res_rev)

## p=0.658 (Not different)


## Run MCMC with rev constraint applied and set XA = ZA
## fit_bay(f_lik_xzeq, f_sp_chr, "output/results/chromosome-fish.rds")

## Read in fusion chain
f_mcmc_chr <- readRDS("output/results/chromosome-fish.rds")
## Remove burnin
f_mcmc_chr <- f_mcmc_chr[-seq_len(10000),]

## Differences between YA and XA
f_mcmc_chr$xy <- f_mcmc_chr$q13 - f_mcmc_chr$q14

fig_fusion_rates_xy <- function(x, taxa){
    profiles.plot(x["xy"], col.line=col.taxa[taxa], opacity=0.9,
                  xlab="Difference between YA and XA/ZA fusion rates",
                  frame.plot=FALSE)
    axis(side=1)
    abline(v=0, lwd=2, lty=2, col=col.vlin)
}

fig_fusion_rates_xy(f_mcmc_chr, "fish")

## What proportion of posterior supports YA > XA/ZA
length(which(f_mcmc_chr$xy > 0))/nrow(f_mcmc_chr)
## 99.5% of posterior probability supports this conclusion








## ## Squamate chromosome analysis (ML)

## ### Definitions of states 
## 1: XY
## 2: ZW
## 3: YA (XXY)
## 4: ZA (ZWW)
## 5: WA (ZZW)

## Note that XYY is not present

## Baseline constraints (biologically enforced)

## 1. Disallow transtions from all fused to other fused states
## 2. Disallow transitions from XY <-> ZW fused / ZW <-> XY fused
constrain_chr_squa <- function(lik)
    constrain(lik, q34~0, q35~0, q43~0, q45~0, q53~0, q54~0,
              q15~0, q14~0, q51~0, q41~0, q23~0, q32~0)

## 3. Constrain rev transtions to all have equal rates
constrain_rev_squa <- function(lik)
    constrain(lik, q42~q31, q52~q31)

## 4a. Constrain YA fusions to be equal to WA
constrain_xyeq_squa <- function(lik)
    constrain(lik, q25~q13)

## 4b. Constrain YA fusions to be equal to ZA
constrain_yzeq_squa <- function(lik)
    constrain(lik, q24~q13)

## 4c. Constrain ZA fusions to be equal to XA
constrain_xzeq_squa <- function(lik)
    constrain(lik, q25~q24)


## define starting points
s_sp_chr <- starting.point.musse(s_chr$phy, k=5)

## build likelihood function
s_lik_chr <- make.mkn(s_chr$phy, s_chr$data, k=5)

## Baseline -- N.B. This is not possible to estimate with ML
s_lik_chr <- constrain_chr_squa(s_lik_chr)
## s_res_chr <- fit_mle(s_lik_chr, s_sp_chr)

## Rev transitions
s_lik_rev <- constrain_rev_squa(s_lik_chr)
s_res_rev <- fit_mle(s_lik_rev, s_sp_chr)

## anova(s_res_rev, s_res_chr)

## Assume that constraining reverse transtions is ok

## YA different from XA?
s_lik_xyeq <- constrain_xyeq_squa(s_lik_rev)
s_res_xyeq <- fit_mle(s_lik_xyeq, s_sp_chr)

anova(s_res_xyeq, s_res_rev)

## p=2.94e-5**

## YA different from ZA?
s_lik_yzeq <- constrain_yzeq_squa(s_lik_rev)
s_res_yzeq <- fit_mle(s_lik_yzeq, s_sp_chr)

anova(s_res_yzeq, s_res_rev)

## p=2.94e-5**

## XA different from ZA?
s_lik_xzeq <- constrain_xzeq_squa(s_lik_rev)
s_res_xzeq <- fit_mle(s_lik_xzeq, s_sp_chr)

anova(s_res_xzeq, s_res_rev)

## likelihoods identical. Ignore

## Run MCMC with only rev constraints applied and setting XA = ZA
## fit_bay(s_lik_xzeq, s_sp_chr, "output/results/chromosome-squa.rds")


## Read in fusion chain
s_mcmc_chr <- readRDS("output/results/chromosome-squa.rds")
## Remove burnin
s_mcmc_chr <- s_mcmc_chr[-seq_len(10000),]

## Differences between YA and XA
s_mcmc_chr$xy <- s_mcmc_chr$q13 - s_mcmc_chr$q24

fig_fusion_rates_xy(s_mcmc_chr, "squa")

## What proportion of posterior supports YA > XA/ZA
length(which(s_mcmc_chr$xy > 0))/nrow(s_mcmc_chr)
## 99.9% of posterior probability supports this conclusion

if (!interactive()){
    to_pdf("output/figs/phylo-fish.pdf", 8, 8,
           fig_phylogeny(f_dat))

    to_pdf("output/figs/phylo-squa.pdf", 8, 8,
           fig_phylogeny(s_dat))

    to_pdf("output/figs/karyotype-fusion-fish.pdf", 4, 4,
           fig_fusion_rates(f_mcmc, "fish"))

    to_pdf("output/figs/karyotype-fusion-squa.pdf", 4, 4,
           fig_fusion_rates(s_mcmc, "squa"))

    to_pdf("output/figs/karyotype-fusion-squa-6par.pdf", 4, 4,
           fig_fusion_rates(s_mcmc_rev, "squa"))

    to_pdf("output/figs/karyotype-residency-squa-6par.pdf", 4.75, 4.75,
           fig_fusion_resid(s_mcmc_rev, "squa"))

    to_pdf("output/figs/chromosome-fusion-fish.pdf", 4, 4,
           fig_fusion_rates_xy(f_mcmc_chr, "fish"))

    to_pdf("output/figs/chromosome-fusion-squa.pdf", 4, 4,
           fig_fusion_rates_xy(s_mcmc_chr, "squa"))

}
    









         



