## function for reading in data and processing it for XY/ZW analysis
process_karyotype <- function(td){
    phy <- td$phy
    dat <- data.frame(td$data)

    ## reformat data for using regular MK models
    states <- sapply(seq_len(nrow(dat)), function(x) paste(dat[x,1], dat[x,2], sep=""))
    states <- factor(states, levels=c("00","01","10","11"), labels=c(1:4))
    states <- as.numeric(states)
    names(states) <- rownames(dat)
    list(phy=phy, states=states)
}


fish_karyotype_dat <- function()
    readRDS("output/data/fish.rds")

squa_karyotype_dat <- function()
    readRDS("output/data/squa.rds")

fish_chromosome_dat <- function()
    readRDS("output/data/fish-indchrom.rds")

squa_chromosome_dat <- function()
    readRDS("output/data/squa-indchrom.rds")

## function for getting diversitree cached objected
cache.mcmc <- function(x)
    diversitree:::get.cache(get.likelihood(x))


## function for loading in an processing raw data
raw_karyotype_dat <- function(){
    dat.all <- read.csv("datasets/vert.data-may19.csv")

    ## fix species names
    dat.all$species <- sapply(dat.all$species, function(x)
                        gsub("[-]", replacement="", x))
    dat.all$species <- sapply(dat.all$species, function(x)
                        gsub("*([0-9])", replacement="", x))

    ## add new column with full name
    dat.all$binom <- sapply(seq_len(nrow(dat.all)), function(x) paste(dat.all[x,"Genus"], dat.all[x, "species"], sep="_"))

    ## change notes column heading
    n <- grep("notes", colnames(dat.all))[1]
    colnames(dat.all)[n] <- "notes"

    ## change karyotype column heading
    k <- grep("Karyotype", colnames(dat.all))[1]
    colnames(dat.all)[k] <- "karyotype"

    dat.all
}


to_pdf <- function(filename, width, height, expr,
                   ..., pointsize=12, verbose=TRUE) {
  if (verbose)
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, width=width, height=height, pointsize=pointsize, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}
