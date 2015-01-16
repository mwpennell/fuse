## generate dataset for hermaphrodite-gonochoristic transition analysis -- fish

## source files for pruning trees
source("R/phyutility.R")

## read in dataset
d.all <- read.csv("datasets/vert.data-may19.csv")

d <- subset(d.all, Higher.taxonomic.group == "Mammalia")

## fix species names
d$species <- sapply(d$species, function(x) gsub("[-]", replacement="", x))
d$species <- sapply(d$species, function(x) gsub("*([0-9])", replacement="", x))

## add new column with full name
d$name <- sapply(seq_len(nrow(d)), function(x) paste(d[x,"Genus"], d[x, "species"], sep="_"))

k <- grep("Karyotype", colnames(d))[1]
n <- grep("notes", colnames(d))[1]
h.names <- c("homomorphic", "XY", "ZW", "complex XY", "complex ZW")
tmp <- which(d[,k] %in% h.names)
dd <- d[tmp,]

states <- rep(1, nrow(dd))
ya <- grep("X2", dd[,n])
states[ya] <- 2
xa <- grep("Y2", dd[,n])
states[xa] <- 3


## build dummy matrix for analysis
hg <- cbind(states, rep(0, length(states)))
rownames(hg) <- dd$name
colnames(hg) <- c("F", "none")
    
## get rid of duplicate rownames
hg <- prune.doubles(hg)

## remove species with NA for one of the traits
hg <- hg[!is.na(hg[,1]),]
hg <- hg[!is.na(hg[,2]),]

## read in tree
t <- read.tree("datasets/mamm.tre")

## match tips to data with sampling algorithm
td <- match.tips.genus(t,hg)

td$data <- td$data[,1]
   

## write RDS file
saveRDS(td, file="output/data/mamm.rds")
