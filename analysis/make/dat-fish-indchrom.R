## generate dataset for hermaphrodite-gonochoristic transition analysis -- fish

## source files for pruning trees
source("R/phyutility.R")

## read in dataset
d.all <- read.csv("datasets/vert.data-may19.csv")

d <- subset(d.all, Higher.taxonomic.group == "Fish")

## fix species names
d$species <- sapply(d$species, function(x) gsub("[-]", replacement="", x))
d$species <- sapply(d$species, function(x) gsub("*([0-9])", replacement="", x))

## add new column with full name
d$name <- sapply(seq_len(nrow(d)), function(x) paste(d[x,"Genus"], d[x, "species"], sep="_"))

## fix column name for notes column
notes <- grep("notes", colnames(d))
colnames(d)[notes] <- "notes"

k <- grep("Karyotype", colnames(d))[1]
h.names <- c("homomorphic", "XY", "ZW", "complex XY", "complex ZW")
tmp <- which(d[,k] %in% h.names)
dd <- d[tmp,]

## XY homomorph or heteromorphs
g <- grep("Genotypic", colnames(dd))[1]
states <- factor(dd[,g])
g.names <- c("male heterogametic", "female heterogametic")
states <- factor(states, levels=g.names, labels=c(1,2))
states <- as.numeric(levels(states)[as.integer(states)])

## Y-A fusions: 3
states[grep("X2", dd$notes)] <- 3

## X-A fusions: 4
states[grep("Y2", dd$notes)] <- 4 

## Z-A fusions: 5
states[grep("W2", dd$notes)] <- 5

## W-A fusions: 6
states[grep("Z2", dd$notes)] <- 6

tmp <- cbind(states, rep(0, length(states)))
rownames(tmp) <- dd$name

tmp <- prune.doubles(tmp)

tmp <- tmp[!is.na(tmp[,1]),]

t <- read.tree("datasets/fish.tre")

td <- match.tips.genus(t,tmp)

td$data <- td$data[,1]

## write RDS file
saveRDS(td, file="output/data/fish-indchrom.rds")







