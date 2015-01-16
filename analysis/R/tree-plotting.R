## Code for plotting fusions

library(diversitree)

col.tree <- list(G=c("#86B8B1","#FA2A00"), K=c("white", "#3D1C00","#F2D694"))


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

notes <- grep("notes", colnames(dat.all))
colnames(dat.all)[notes] <- "notes"



## ## Fish
dat.fish <- subset(dat.all, Higher.taxonomic.group == "Fish")

## How many Y-A fusions
ya.fish <- grep("X2", dat.fish$notes)
names(ya.fish) <- dat.fish$binom[ya.fish]
het.y.fish <- sapply(names(ya.fish), function(x) {strsplit(x, split=" ")[[1]][1]})

## How many X-A fusions
xa.fish <- grep("Y2", dat.fish$notes)
names(xa.fish) <- dat.fish$binom[xa.fish]
hom.x.fish <- sapply(names(xa.fish), function(x) {strsplit(x, split=" ")[[1]][1]})

## How many Z-A fusions
za.fish <- grep("W2", dat.fish$notes)
names(za.fish) <- dat.fish$binom[za.fish]
hom.z.fish <- sapply(names(za.fish), function(x) {strsplit(x, split=" ")[[1]][1]})

## How many W-A fusions
wa.fish <- grep("Z2", dat.fish$notes)
names(wa.fish) <- dat.fish$binom[wa.fish]
het.w.fish <- sapply(names(wa.fish), function(x) {strsplit(x, split=" ")[[1]][1]})

## read in tree matched data
fdat <- readRDS("output/data/fish.rds")

fphy <- fdat$phy
fdat <- as.data.frame(fdat$data)

for (i in seq_len(nrow(fdat))){
    if (fdat[i,"K"] == 1){
        if (rownames(fdat)[i] %in% c(hom.x.fish, hom.z.fish))
            fdat[i, "K"] <- 2
    }
}

genus.fish <- sub("_.+$", "", fphy$tip.label)
trait.plot(ladderize(fphy), dat=fdat, cols=col.tree,
           lab=c("Genotype", "Fused"),
           str=list(c("XY", "ZW"), c("None", "Het. Chrom", "Hom. Chrom")),
           legend=FALSE,
           class=genus.fish, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)


## ## Squamates
dat.squa <- subset(dat.all, Order == "Squamata")

## How many Y-A fusions
ya.squa <- grep("X2", dat.squa$notes)
names(ya.squa) <- dat.squa$binom[ya.squa]
het.y.squa <- sapply(names(ya.squa), function(x) {strsplit(x, split=" ")[[1]][1]})

## How many X-A fusions
xa.squa <- grep("Y2", dat.squa$notes)
names(xa.squa) <- dat.squa$binom[xa.squa]
hom.x.squa <- sapply(names(xa.squa), function(x) {strsplit(x, split=" ")[[1]][1]})

## How many Z-A fusions
za.squa <- grep("W2", dat.squa$notes)
names(za.squa) <- dat.squa$binom[za.squa]
hom.z.squa <- sapply(names(za.squa), function(x) {strsplit(x, split=" ")[[1]][1]})

## How many W-A fusions
wa.squa <- grep("Z2", dat.squa$notes)
names(wa.squa) <- dat.squa$binom[wa.squa]
het.w.squa <- sapply(names(wa.squa), function(x) {strsplit(x, split=" ")[[1]][1]})

## read in tree matched data
sdat <- readRDS("output/data/squa.rds")

sphy <- sdat$phy
sdat <- as.data.frame(sdat$data)

for (i in seq_len(nrow(sdat))){
    if (sdat[i,"K"] == 1){
        if (rownames(sdat)[i] %in% c(hom.x.squa, hom.z.squa))
            sdat[i, "K"] <- 2
    }
}

genus.squa <- sub("_.+$", "", sphy$tip.label)
trait.plot(ladderize(sphy), dat=sdat, cols=col.tree,
           lab=c("Genotype", "Fused"),
           str=list(c("XY", "ZW"), c("None", "Het. Chrom", "Hom. Chrom")),
           legend=FALSE,
           class=genus.squa, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)


dev.off()

pdf("output/figs/tree-fig.pdf", width=16, height=8)
par(mfrow=c(1,2))
par(oma=c(0.75,0.75,0.75,0.75))
trait.plot(ladderize(fphy), dat=fdat, cols=col.tree,
           lab=c("Genotype", "Fused"),
           str=list(c("XY", "ZW"), c("None", "Het. Chrom", "Hom. Chrom")),
           legend=FALSE,
           class=genus.fish, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)
trait.plot(ladderize(sphy), dat=sdat, cols=col.tree,
           lab=c("Genotype", "Fused"),
           str=list(c("XY", "ZW"), c("None", "Het. Chrom", "Hom. Chrom")),
           legend=FALSE,
           class=genus.squa, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)
dev.off()



