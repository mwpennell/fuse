#!/usr/bin/env Rscript

pkgs.cran <- c(knitr="1.5")

pkgs.github <- c("mwpennell/arbutus"="1.0",
                 "richfitz/sowsear"="0.1-1")

pkgs <- installed.packages()

to.install <- setdiff(names(pkgs.cran), rownames(pkgs))

installed <- intersect(names(pkgs.cran), rownames(pkgs))
to.upgrade <-
  installed[numeric_version(pkgs[installed, "Version"]) <
            numeric_version(pkgs.cran[installed])]

msg.install <- (if (length(to.install) == 0) character(0) else
                paste("installing:", paste(to.install, collapse=", ")))
msg.upgrade <- (if (length(to.upgrade) == 0) character(0) else
                paste("upgrading:", paste(to.upgrade, collapse=", ")))
msg <- c(msg.install, msg.upgrade)
if (length(msg) > 0) {
  message(paste(msg, collapse="\n"))
  install.packages(c(to.install, to.upgrade))
}

# GitHub:
pkgs.github.short <- sub("^.+/", "", names(pkgs.github))

tr <- structure(names(pkgs.github), names=pkgs.github.short)

to.install <- setdiff(pkgs.github.short, rownames(pkgs))
installed <- intersect(pkgs.github.short, rownames(pkgs))

to.upgrade <-
  installed[numeric_version(pkgs[installed, "Version"]) <
            numeric_version(pkgs.github[tr[installed]])]

msg.install <- (if (length(to.install) == 0) character(0) else
                paste("installing:", paste(to.install, collapse=", ")))
msg.upgrade <- (if (length(to.upgrade) == 0) character(0) else
                paste("upgrading:", paste(to.upgrade, collapse=", ")))
msg <- c(msg.install, msg.upgrade)
if (length(msg) > 0) {
  message(paste(msg, collapse="\n"))
  install.packages("devtools")
  library(devtools)
  for (pkg in tr[c(to.install, to.upgrade)]) {
    install_github(pkg)
  }
}



