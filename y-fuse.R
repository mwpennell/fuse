## # Analysis of sex chromosome fusions

## ### Preliminaries

### Set options
##+ echo=FALSE, results=FALSE
knitr::opts_chunk$set(tidy=FALSE)

## load diversitree
library(diversitree)

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

## function for getting diversitree cached objected
cache.mcmc <- function(x)
    diversitree:::get.cache(get.likelihood(x))

## set colors
col.tree <- list(G=c("#86B8B1","#FA2A00"), K=c("#F2D694","#3D1C00"))

## ### Plot sex chromosomes on tree
## Fish data
fdat <- readRDS("output/data/fish.rds")
genus.fish <- sub("_.+$", "", fdat$phy$tip.label)
trait.plot(ladderize(fdat$phy), dat=as.data.frame(fdat$dat), cols=col.tree,
           lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")),
           class=genus.fish, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)

## Number of species in analysis
nrow(fdat$dat)


## Squamate data 
sdat <- readRDS("output/data/squa.rds")
genus.squa <- sub("_.+$", "", sdat$phy$tip.label)
trait.plot(ladderize(sdat$phy), dat=as.data.frame(sdat$dat), cols=col.tree,
           lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")),
           class=genus.squa, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)


## Number of species in analysis
nrow(sdat$dat)


## # Part A: Model fitting

## ## Model 1:

## 4 parameters:
## 1. transition rate to XY from any other state
## 2. transition rate to ZW from any other state
## 3. rate of fusion in XY lineages
## 4. rate of fusion in ZW lineages

## Read in data

feq <- readRDS("output/results/fish-4par.rds")
seq <- readRDS("output/results/squa-4par.rds")

## Get cache for likelihood

feq.cache <- cache.mcmc(feq)
seq.cache <- cache.mcmc(seq)

## Remove burning
feq <- feq[-seq_len(10000),]
seq <- seq[-seq_len(10000),]

## Are rates of fusion different between XY and ZW

## ### Fish
feq$ff <- sapply(seq_len(nrow(feq)), function(x)
                 {feq[x,"q12"] - feq[x,"q34"]})

profiles.plot(feq["ff"], col.line="#053061")

## posterior density greater than 0
length(which(feq$ff >= 0))/length(feq$ff)



## ### Squamates
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


falt <- readRDS("output/results/fish-5par.rds")
salt <- readRDS("output/results/squa-5par.rds")

## Get cache for likelihood

falt.cache <- cache.mcmc(falt)
salt.cache <- cache.mcmc(salt)

## Remove burning
falt <- falt[-seq_len(10000),]
salt <- salt[-seq_len(10000),]

## Are rates of fusion different between XY and ZW

## ### Fish

falt$ff <- sapply(seq_len(nrow(falt)), function(x)
                 {falt[x,"q12"] - falt[x,"q34"]})

profiles.plot(falt["ff"], col.line="#053061")

## posterior density greater than 0
length(which(falt$ff >= 0))/length(falt$ff)



## ### Squamates

salt$ff <- sapply(seq_len(nrow(salt)), function(x)
                 {salt[x,"q12"] - salt[x,"q34"]})

profiles.plot(salt["ff"], col.line="#053061")

length(which(salt$ff >= 0))/length(salt$ff)


## ## Model 3
## 6 parameters
## 1. XY fusion rate
## 2. ZX fusion rate
## 3. XY -> ZW
## 4. ZX -> XY
## 5. XXY -> XY
## 6. XXY -> ZW
fsix <- readRDS("output/results/fish-6par.rds")
ssix <- readRDS("output/results/squa-6par.rds")

## Get cache for likelihood

fsix.cache <- cache.mcmc(fsix)
ssix.cache <- cache.mcmc(ssix)

## Remove burning
fsix <- fsix[-seq_len(10000),]
ssix <- ssix[-seq_len(10000),]

## Are rates of fusion different between XY and ZW

## ### Fish

fsix$ff <- sapply(seq_len(nrow(fsix)), function(x)
                 {fsix[x,"q12"] - fsix[x,"q34"]})

profiles.plot(fsix["ff"], col.line="#053061")

## posterior density greater than 0
length(which(fsix$ff >= 0))/length(fsix$ff)



## ### Squamates

ssix$ff <- sapply(seq_len(nrow(ssix)), function(x)
                 {ssix[x,"q12"] - ssix[x,"q34"]})

profiles.plot(ssix["ff"], col.line="#053061")

length(which(ssix$ff >= 0))/length(ssix$ff)


six <- data.frame(fish=fsix$ff, squa=ssix$ff)
par(mar=c(5, 4, 2, 2))
profiles.plot(six[c("fish", "squa")], col.line=c("#86B8B1","#FA2A00"),
              xlim=c(-0.01, 0.05),
              xlab="Difference in fusion rates between XY and ZW",
              opacity=0.7, yaxt="n", frame.plot=FALSE, cex=2)
abline(v=0, lwd=2, lty=2, col="#3D1C00")




## # Part B: Counting fusions
## Read in raw datastes
dat.all <- read.csv("datasets/vert.data-may19.csv")

## fix species names
dat.all$species <- sapply(dat.all$species, function(x)
                        gsub("[-]", replacement="", x))
dat.all$species <- sapply(dat.all$species, function(x)
                        gsub("*([0-9])", replacement="", x))

## add new column with full name
dat.all$binom <- sapply(seq_len(nrow(dat.all)), function(x) paste(dat.all[x,"Genus"], dat.all[x, "species"], sep="_"))

## change karyotype column heading
k <- grep("Karyotype", colnames(dat.all))[1]
colnames(dat.all)[k] <- "karyotype"



## ## Fish
dat.fish <- subset(dat.all, Higher.taxonomic.group == "Fish")

## How many Y-A fusions
ya.fish <- grep("X2", dat.fish$notes)
names(ya.fish) <- dat.fish$binom[ya.fish]
tmp <- sapply(names(ya.fish), function(x) {strsplit(x, split=" ")[[1]][1]})
n.ya.fish <- length(unique(tmp))
n.ya.fish ## 42

## How many X-A fusions
xa.fish <- grep("Y2", dat.fish$notes)
names(xa.fish) <- dat.fish$binom[xa.fish]
tmp <- sapply(names(xa.fish), function(x) {strsplit(x, split=" ")[[1]][1]})
n.xa.fish <- length(unique(tmp))
n.xa.fish ## 3

## How many Z-A fusions
za.fish <- grep("W2", dat.fish$notes)
names(za.fish) <- dat.fish$binom[za.fish]
tmp <- sapply(names(za.fish), function(x) {strsplit(x, split=" ")[[1]][1]})
n.za.fish <- length(unique(tmp))
n.za.fish ## 3

## How many W-A fusions
wa.fish <- grep("Z2", dat.fish$notes)
names(wa.fish) <- dat.fish$binom[wa.fish]
tmp <- sapply(names(wa.fish), function(x) {strsplit(x, split=" ")[[1]][1]})
n.wa.fish <- length(unique(tmp))
n.wa.fish ## 1

## total number of species with karyotype information
n.tot.fish <- length(which(dat.fish$karyotype != ""))
n.tot.fish

## unique species only
n.uni.fish <- length(unique(dat.fish[which(dat.fish$karyotype != ""),"binom"]))
n.uni.fish

## number of families
n.fam.fish <- length(unique(dat.fish[which(dat.fish$karyotype != ""),"Family"]))
n.fam.fish

## ## Squamates
dat.squa <- subset(dat.all, Order=="Squamata")

## How many Y-A fusions
ya.squa <- grep("X2", dat.squa$notes)
names(ya.squa) <- dat.squa$binom[ya.squa]
tmp <- sapply(names(ya.squa), function(x) {strsplit(x, split=" ")[[1]][1]})
n.ya.squa <- length(unique(tmp))
n.ya.squa ## 40

## How many X-A fusions
xa.squa <- grep("Y2", dat.squa$notes)
names(xa.squa) <- dat.squa$binom[xa.squa]
tmp <- sapply(names(xa.squa), function(x) {strsplit(x, split=" ")[[1]][1]})
n.xa.squa <- length(unique(tmp))
n.xa.squa ## 0

## How many Z-A fusions
za.squa <- grep("W2", dat.squa$notes)
names(za.squa) <- dat.squa$binom[za.squa]
tmp <- sapply(names(za.squa), function(x) {strsplit(x, split=" ")[[1]][1]})
n.za.squa <- length(unique(tmp))
n.za.squa ## 4

## How many W-A fusions
wa.squa <- grep("Z2", dat.squa$notes)
names(wa.squa) <- dat.squa$binom[wa.squa]
tmp <- sapply(names(wa.squa), function(x) {strsplit(x, split=" ")[[1]][1]})
n.wa.squa <- length(unique(tmp))
n.wa.squa ## 0

## total number of species with karyotype information
n.tot.squa <- length(which(dat.squa$karyotype != ""))
n.tot.squa

## unique species only
n.uni.squa <- length(unique(dat.squa[which(dat.squa$karyotype != ""),"binom"]))
n.uni.squa

## number of families
n.fam.squa <- length(unique(dat.squa[which(dat.squa$karyotype != ""),"Family"]))
n.fam.squa



## ## Turtles

dat.chel <- subset(dat.all, Order=="Chelonia")

## How many Y-A fusions
ya.chel <- grep("X2", dat.chel$notes)
names(ya.chel) <- dat.chel$binom[ya.chel]
n.ya.chel <- length(ya.chel)
n.ya.chel ## 0

## How many X-A fusions
xa.chel <- grep("Y2", dat.chel$notes)
names(xa.chel) <- dat.chel$binom[xa.chel]
n.xa.chel <- length(xa.chel)
n.xa.chel ## 0

## How many Z-A fusions
za.chel <- grep("W2", dat.chel$notes)
names(za.chel) <- dat.chel$binom[za.chel]
n.za.chel <- length(za.chel)
n.za.chel ## 0

## How many W-A fusions
wa.chel <- grep("Z2", dat.chel$notes)
names(wa.chel) <- dat.chel$binom[wa.chel]
n.wa.chel <- length(wa.chel)
n.wa.chel ## 0

## total number of species with karyotype information
n.tot.chel <- length(which(dat.chel$karyotype != ""))
n.tot.chel

## unique species only
n.uni.chel <- length(unique(dat.chel[which(dat.chel$karyotype != ""),"binom"]))
n.uni.chel



## ## Mammals

dat.mamm <- subset(dat.all, Higher.taxonomic.group=="Mammalia")

## How many Y-A fusions
ya.mamm <- grep("X2", dat.mamm$notes)
names(ya.mamm) <- dat.mamm$binom[ya.mamm]
tmp <- sapply(names(ya.mamm), function(x) {strsplit(x, split=" ")[[1]][1]})
n.ya.mamm <- length(unique(tmp))
n.ya.mamm ## 20

## How many X-A fusions
xa.mamm <- grep("Y2", dat.mamm$notes)
names(xa.mamm) <- dat.mamm$binom[xa.mamm]
tmp <- sapply(names(xa.mamm), function(x) {strsplit(x, split=" ")[[1]][1]})
n.xa.mamm <- length(unique(tmp))
n.xa.mamm ## 25

## How many Z-A fusions
za.mamm <- grep("W2", dat.mamm$notes)
names(za.mamm) <- dat.mamm$binom[za.mamm]
n.za.mamm <- length(za.mamm)
n.za.mamm ## 0

## How many W-A fusions
wa.mamm <- grep("Z2", dat.mamm$notes)
names(wa.mamm) <- dat.mamm$binom[wa.mamm]
n.wa.mamm <- length(wa.mamm)
n.wa.mamm ## 0

## total number of species with karyotype information
n.tot.mamm <- length(which(dat.mamm$karyotype != ""))
n.tot.mamm

## unique species only
n.uni.mamm <- length(unique(dat.mamm[which(dat.mamm$karyotype != ""),"binom"]))
n.uni.mamm



## ## Birds

dat.aves <- subset(dat.all, Higher.taxonomic.group=="Aves")

## How many Y-A fusions
ya.aves <- grep("X2", dat.aves$notes)
names(ya.aves) <- dat.aves$binom[ya.aves]
n.ya.aves <- length(ya.aves)
n.ya.aves ## 0

## How many X-A fusions
xa.aves <- grep("Y2", dat.aves$notes)
names(xa.aves) <- dat.aves$binom[xa.aves]
n.xa.aves <- length(xa.aves)
n.xa.aves ## 0

## How many Z-A fusions
za.aves <- grep("W2", dat.aves$notes)
names(za.aves) <- dat.aves$binom[za.aves]
n.za.aves <- length(za.aves)
n.za.aves ## 0

## How many W-A fusions
wa.aves <- grep("Z2", dat.aves$notes)
names(wa.aves) <- dat.aves$binom[wa.aves]
n.wa.aves <- length(wa.aves)
n.wa.aves ## 0

## total number of species with karyotype information
n.tot.aves <- length(which(dat.aves$karyotype != ""))
n.tot.aves

## unique species only
n.uni.aves <- length(unique(dat.aves[which(dat.aves$karyotype != ""),"binom"]))
n.uni.aves


## ## Amphibians

dat.amph <- subset(dat.all, Higher.taxonomic.group=="Amphibia")

## How many Y-A fusions
ya.amph <- grep("X2", dat.amph$notes)
names(ya.amph) <- dat.amph$binom[ya.amph]
n.ya.amph <- length(ya.amph)
n.ya.amph ## 0

## How many X-A fusions
xa.amph <- grep("Y2", dat.amph$notes)
names(xa.amph) <- dat.amph$binom[xa.amph]
n.xa.amph <- length(xa.amph)
n.xa.amph ## 0

## How many Z-A fusions
za.amph <- grep("W2", dat.amph$notes)
names(za.amph) <- dat.amph$binom[za.amph]
n.za.amph <- length(za.amph)
n.za.amph ## 0

## How many W-A fusions
wa.amph <- grep("Z2", dat.amph$notes)
names(wa.amph) <- dat.amph$binom[wa.amph]
n.wa.amph <- length(wa.amph)
n.wa.amph ## 0

## total number of species with karyotype information
n.tot.amph <- length(which(dat.amph$karyotype != ""))
n.tot.amph

## unique species only
n.uni.amph <- length(unique(dat.amph[which(dat.amph$karyotype != ""),"binom"]))
n.uni.amph


## ## Overall results
ya <- c(n.ya.fish, n.ya.squa, n.ya.chel, n.ya.mamm, n.ya.aves, n.ya.amph)
xa <- c(n.xa.fish, n.xa.squa, n.xa.chel, n.xa.mamm, n.xa.aves, n.xa.amph)
wa <- c(n.wa.fish, n.wa.squa, n.wa.chel, n.wa.mamm, n.wa.aves, n.wa.amph)
za <- c(n.za.fish, n.za.squa, n.za.chel, n.za.mamm, n.za.aves, n.za.amph)
tot <- c(n.tot.fish, n.tot.squa, n.tot.chel, n.tot.mamm, n.tot.aves, n.tot.amph)
uni <- c(n.uni.fish, n.uni.squa, n.uni.chel, n.uni.mamm, n.uni.aves, n.uni.amph)

fusion <- data.frame(tot, uni, ya, xa, wa, za)
colnames(fusion) <- c("Tot.records", "Unique.spp", "YA", "XA", "WA", "ZA") 
rownames(fusion) <- c("Fish", "Squamates", "Turtles", "Mammals",
                      "Birds", "Amphibians")

fusion




if (!interactive()){
    dev.off()
    pdf("output/figs/trait-phylo-fish.pdf", width=8, height=8)
    par(oma=c(0.75,0.75,0.75,0.75))
    trait.plot(ladderize(fdat$phy), dat=as.data.frame(fdat$dat),
               cols=col.tree, lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"),
               c("No", "Yes")), class=genus.fish, quiet=TRUE, w=1/10,
               margin=1/1.4, cex.legend = 1)
    dev.off()

    pdf("output/figs/trait-phylo-squamates.pdf", width=8, height=8)
    par(oma=c(0.75,0.75,0.75,0.75))
    trait.plot(ladderize(sdat$phy), dat=as.data.frame(sdat$dat),
               cols=col.tree, lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"),
               c("No", "Yes")), class=genus.squa, quiet=TRUE, w=1/10,
               margin=1/1.4, cex.legend = 1)
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

    pdf("output/figs/fusion-xyzw.pdf")
    par(mar=c(5, 4, 2, 2))
    profiles.plot(six[c("fish", "squa")], col.line=c("#86B8B1","#FA2A00"),
               xlim=c(-0.01, 0.05),
               xlab="Difference in fusion rates between XY and ZW",
               opacity=0.7, yaxt="n", frame.plot=FALSE, cex=2)
    abline(v=0, lwd=2, lty=2, col="#3D1C00")
    dev.off()

    write.csv(fusion, "output/results/fusion-table.csv")
    
  
}

