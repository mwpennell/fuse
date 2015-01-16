library(diversitree)

## set colors for the two characters
col.tree <- list(G=c("#86B8B1","#FA2A00"), K=c("#F2D694","#3D1C00"))

## ### Plot sex chromosomes on tree
## Fish data
fdat <- readRDS("output/data/fish.rds")

phy <- fdat$phy
dat <- as.data.frame(fdat$data)

## Get genera names
genus.fish <- sub("_.+$", "", phy$tip.label)

## Plot the tree
trait.plot(ladderize(phy), dat=dat, cols=col.tree,
           lab=c("Genotype", "Fused"), str=list(c("XY", "ZW"), c("No", "Yes")),
           class=genus.fish, quiet=TRUE, w=1/10, margin=1/1.4, cex.legend = 1)
