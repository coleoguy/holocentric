library(phytools)
library(chromePlus)
trees <- read.nexus("misof.backbone.nex")
trees <- read.nexus("rainford.backbone.nex")
dat <- read.csv("data.invert.csv", as.is=T)


source("functions.R")

foo <- getData2(trees, dat)
tree <- foo[[1]][[1]]
dat <- foo[[2]]
# mono is zero

chrom.type <- dat$chrom
names(chrom.type) <- dat$genus

maps.eq <- make.simmap(tree = tree,
                    x = chrom.type,
                    model = "ARD",
                    pi = "equal",
                    nsim=100)
maps.est <- make.simmap(tree = tree,
                    x = chrom.type,
                    model = "ARD",
                    pi = "estimated",
                    nsim=100)
par(mfcol=c(1,2))
smap.eq <- densityMap(trees = maps.eq,
                   res=50, lwd=.75,
                   check=FALSE,
                   type="phylogram",
                   direction="rightwards")
smap.est <- densityMap(trees = maps.est,
                   res=50, lwd=.75,
                   check=FALSE,
                   type="fan",
                   direction="rightwards",
                   fsize=.25)


