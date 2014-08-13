## Analysis of sex chromosome fusions
library(diversitree)

## Some helper functions

## function for getting 95% HPD from mcmc.samples
hdr <- diversitree:::hdr

## function for determining what proportion of 95% HPD overlaps with 0
z.hdr <- function(x){
   hpd <- hdr(x)
   ## all samples in 95% HPD
   nf <- which(which(x >= hpd.lim[1]) %in% which(x <= hpd.lim[2]))
   ## lower than 0
   l <- length(which(x[nf] < 0))/length(x)
   ## greater than 0
   g <- length(which(x[nf] > 0))/length(x)
   list(lower.than.zero=l, greater.than.zero=g)
}

cache.mcmc <- function(x)
    diversitree:::get.cache(get.likelihood(x))

## set colors
col.tree <- list(G=c("#86B8B1","#FA2A00"), K=c("#F2D694","#3D1C00"))

## First, look at actual data
fdat <- readRDS("output/data/fish.rds")
genus.fish <- sub("_.+$", "", fdat$phy$tip.label)
trait.plot(ladderize(fdat$phy), dat=as.data.frame(fdat$dat), cols=col.tree, lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")), class=genus.fish, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)

sdat <- readRDS("output/data/squa.rds")
genus.squa <- sub("_.+$", "", sdat$phy$tip.label)
trait.plot(ladderize(sdat$phy), dat=as.data.frame(sdat$dat), cols=col.tree, lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")), class=genus.squa, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)


## ## Model 1:

## 4 parameters:
## 1. transition rate to XY from any other state
## 2. transition rate to ZW from any other state
## 3. rate of fusion in XY lineages
## 4. rate of fusion in ZW lineages

## Read in data

feq <- readRDS("output/results/fish-feq.rds")
seq <- readRDS("output/results/squa-feq.rds")

## Get cache for likelihood

feq.cache <- cache.mcmc(feq)
seq.cache <- cache.mcmc(seq)

## Remove burning
feq <- feq[-seq_len(10000),]
seq <- seq[-seq_len(10000),]

## Are rates of fusion different between XY and ZW

## Fish
feq$fus <- sapply(seq_len(nrow(feq)), function(x)
                  {feq[x,"qxy"] - feq[x,"qzw"]})
profiles.plot(feq["fus"], col.line="#053061")

## posterior density greater than 0
length(which(feq$fus >= 0))/length(feq$fus)


feq$ff <- sapply(seq_len(nrow(feq)), function(x)
                 {feq[x,"q12"] - feq[x,"q34"]})

profiles.plot(feq["ff"], col.line="#053061")

## posterior density greater than 0
length(which(feq$ff >= 0))/length(feq$ff)



## Squamates
seq$fus <- sapply(seq_len(nrow(seq)), function(x)
                  {seq[x,"qxy"] - seq[x,"qzw"]})
profiles.plot(seq["fus"], col.line="#053061")

## posterior density greater than 0
length(which(seq$fus >= 0))/length(seq$fus)


seq$ff <- sapply(seq_len(nrow(seq)), function(x)
                 {seq[x,"q12"] - seq[x,"q34"]})

profiles.plot(seq["ff"], col.line="#053061")

length(which(seq$ff >= 0))/length(seq$ff)




## ## Model 2
## 5 parameters
## 1. XY fusion rate
## 2. ZX fusion rate
## 3. XY -> ZW
## 4. ZX -> XY
## 5. XXY -> XY

## all others set to 0

## 1. disallow XXY to ZW and ZZW to XY. 2. set ZZW to ZW equal to 0
## Read in data

falt <- readRDS("output/results/fish-5par.rds")
salt <- readRDS("output/results/squa-5par.rds")

## Get cache for likelihood

falt.cache <- cache.mcmc(falt)
salt.cache <- cache.mcmc(salt)

## Remove burning
falt <- falt[-seq_len(10000),]
salt <- salt[-seq_len(10000),]

## Are rates of fusion different between XY and ZW

## Fish

falt$ff <- sapply(seq_len(nrow(falt)), function(x)
                 {falt[x,"q12"] - falt[x,"q34"]})

profiles.plot(falt["ff"], col.line="#053061")

## posterior density greater than 0
length(which(falt$ff >= 0))/length(falt$ff)



## Squamates

salt$ff <- sapply(seq_len(nrow(salt)), function(x)
                 {salt[x,"q12"] - salt[x,"q34"]})

profiles.plot(salt["ff"], col.line="#053061")

length(which(salt$ff >= 0))/length(salt$ff)








if (!interactive()){
    dev.off()
    pdf("output/figs/trait-phylo-fish.pdf", width=8, height=8)
    par(oma=c(0.75,0.75,0.75,0.75))
    trait.plot(ladderize(fdat$phy), dat=as.data.frame(fdat$dat), cols=col.tree, lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")), class=genus.fish, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)
    dev.off()

    pdf("output/figs/trait-phylo-squamates.pdf", width=8, height=8)
    par(oma=c(0.75,0.75,0.75,0.75))
    trait.plot(ladderize(sdat$phy), dat=as.data.frame(sdat$dat), cols=col.tree, lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")), class=genus.squa, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)
    dev.off()

    pdf("output/figs/fusion-rates-xy-zw.pdf", width=9, height=5)
    par(mfrow=c(1,2))
    par(mar=c(5,1,4,2))
    profiles.plot(feq["ff"], col.line="#86B8B1", xlab="", yaxt="n",
                  ylab="", opacity = 1)
    abline(v=0, lwd=2, lty=2, col="#3D1C00")
    par(mar=c(5,0,4,2))
    profiles.plot(seq["ff"], col.line="#86B8B1", xlab="", ylab="",
                  yaxt="n",  opacity = 1)
    abline(v=0, lwd=2, lty=2, col="#3D1C00")
    dev.off()
    
  
}

