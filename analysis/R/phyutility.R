library(geiger)
match.tips.genus <- function(phy, data){
    datsp <- rownames(data)
    
    tregn <- sapply(phy$tip.label, function(x) return(strsplit(x, split="_")[[1]][1]))
    datgn <- sapply(datsp, function(x) return(strsplit(x, split="_")[[1]][1]))

    
    ## step 1: drop all tips for which there is no genus match
    todrop <- vector()
    for (i in seq_len(length(tregn))){
        if (!tregn[i] %in% datgn)
            todrop <- c(todrop, phy$tip.label[i])

    }
    prphy <- geiger:::.drop.tip(phy, tip=todrop) ## drop all from tree

    ## update treesp, treegen
    tresp <- prphy$tip.label
    tregn <- tregn[-which(phy$tip.label %in% todrop)]

    ## update phy
    phy <- prphy


    ## step 2: build table of perfect and imperfect matches
    tregn.uni <- unique(tregn)
    mat <- matrix(nrow=length(tregn.uni), ncol=3)
    colnames(mat) <- c("sp.tree", "sp.match", "gn.match")
    rownames(mat) <- tregn.uni

    for (j in seq_len(length(tregn.uni))){

        gn <- as.numeric(which(tregn == tregn.uni[j]))
        mat[j,1] <- length(gn)
        sp.match <- length(which(tresp[gn] %in% datsp))
        mat[j,2] <- sp.match
        gn.match <- length(which(datgn == tregn.uni[j]))
        mat[j,3] <- gn.match
      
    }


    
    ## step 3: find number of possible replacements
    repl <- vector()
    for (k in seq_len(nrow(mat))){
        if (mat[k,2] <= 1 && mat[k,1] > mat[k,2]){
            tmp <- mat[k,3] - mat[k,2]
        } else {
            tmp <- 0
        }
        repl <- c(repl, tmp)
    }
    mat <- cbind(mat, repl)


    
    ## step 4: replace via sampling
    for (m in seq_len(nrow(mat))){
  
        ## if we should replace
        if (mat[m,"repl"] >= 1){
            tips <- tresp[which(tregn == rownames(mat)[m])]
            dat  <- datsp[which(datgn == rownames(mat)[m])]

            
            if (mat[m,"sp.match"] == 1){
                ## eliminate possibility of replacing itself
                tmp <- tips[which(tips %in% dat)]
                tips <- tips[-which(tips == tmp)]
                dat  <- dat[-which(dat == tmp)]

                ## replace one of the remaining tips
                ## choose first one since it does not matter
                phy$tip.label[which(phy$tip.label == tips[1])] <- sample(dat, 1)
                
            } else {
                ## decide if we can replace with 1 or 2
                if (mat[m,"sp.tree"] >= 2 && mat[m, "repl"] >= 2){
                    ## replace 2
                    phy$tip.label[which(phy$tip.label %in% tips[c(1,2)])] <- sample(dat, 2)
                } else {
                    ## replace 1
                    phy$tip.label[which(phy$tip.label == tips[1])] <- sample(dat, 1)
                }
            }
        }
    }

    ## drop remaining tips
    extra <- phy$tip.label[-which(phy$tip.label %in% datsp)]
    newphy <- geiger:::.drop.tip(phy, tip=extra)

    ## run tree data on new tree
    td <- treedata(newphy, data, warnings=FALSE)
}


## function for pruning out duplicates according to rownames
prune.doubles <- function(d){
    dup <- rownames(d)[which(duplicated(rownames(d)))]
    for (i in dup){
        tmp <- which(rownames(d) == i)
        for (j in 1:ncol(d)){
            if (length(unique(d[tmp,j])) != 1) ## if differenent, assign NA
                d[tmp,j] <- rep(NA, length(tmp))
        }
    }
    d[unique(rownames(d)),]
}
