## generate dataset for hermaphrodite-gonochoristic transition analysis -- fish

## source files for pruning trees
source("R/phyutility.R")

## read in dataset
d.all <- read.csv("datasets/vert.data-may19.csv")

d <- subset(d.all, Order == "Squamata")

## fix species names
d$species <- sapply(d$species, function(x) gsub("[-]", replacement="", x))
d$species <- sapply(d$species, function(x) gsub("*([0-9])", replacement="", x))

## add new column with full name
d$name <- sapply(seq_len(nrow(d)), function(x) paste(d[x,"Genus"], d[x, "species"], sep="_"))

k <- grep("Karyotype", colnames(d))[1]
h.names <- c("homomorphic", "XY", "ZW", "complex XY", "complex ZW")
tmp <- which(d[,k] %in% h.names)
dd <- d[tmp,]

k.code <- function(x){
    if (x == "complex XY" | x == "complex ZW"){
        y <- 1
    } else {
        y <- 0
    }
    y
}

k.states <- sapply(dd[,k], function(x) k.code(x))

g <- grep("Genotypic", colnames(dd))[1]
g.states <- factor(dd[,g])
g.names <- c("male heterogametic", "female heterogametic")
g.states <- factor(g.states, levels=g.names, labels=c(0,1))
g.states <- as.numeric(levels(g.states)[as.integer(g.states)])


## build matrix for analysis
hg <- cbind(g.states, k.states)
rownames(hg) <- dd$name
colnames(hg) <- c("G", "K")
    
## get rid of duplicate rownames
hg <- prune.doubles(hg)

## remove species with NA for one of the traits
hg <- hg[!is.na(hg[,1]),]
hg <- hg[!is.na(hg[,2]),]

## read in tree
t <- read.tree("datasets/squa.tre")

## match tips to data with sampling algorithm
td <- match.tips.genus(t,hg)

## write RDS file
saveRDS(td, file="output/data/squa.rds")
