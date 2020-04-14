library(phytools)
library(chromePlus)
library(viridis)
trees <- read.nexus("../data/misof.backbone.nex")
trees <- read.nexus("../data/rainford.backbone.nex")
dat <- read.csv("../data/data.invert.csv", as.is=T)


source("functions.R")

foo <- getData2(trees, dat)
tree <- foo[[1]][[1]]
dat <- foo[[2]]
# mono is zero

chrom.type <- dat$chrom
names(chrom.type) <- dat$genus

maps.est <- make.simmap(tree = tree,
                    x = chrom.type,
                    model = "ARD",
                    pi = "estimated",
                    nsim=100)
smap.est <- densityMap(trees = maps.est,
                   res=50, lwd=.75,
                   check=FALSE,
                   type="fan",
                   direction="rightwards",
                   fsize=.25)
smap.est2 <- setMap(smap.est, c("gray", "black"))
plot(smap.est2, fsize=c(.00001, .4),
     lwd=.95,type="fan")
new.haps <- round(log(dat$haploid)*10)
cols <- viridis(38)
tip.cols <- c()
foo <- c()
for(i in 1:length(tree$tip.label)){
  hit <- new.haps[dat$genus == tree$tip.label[i]]
  foo[i] <- hit
  tip.cols[i] <- cols[(hit-6)]
}
tiplabels(pch=16, col = tip.cols, cex=.45)

library(plotrix)
plot(0,0,col="white", xlim=c(0,10),ylim=c(0,10))
gradient.rect(xleft=2,
              ybottom=2,
              xright=6,
              ytop=4,
              col=cols,
              nslices=50,gradient="x",border=NA)
# export main figure as 6"x6" pdf
# this is the color range
c(2, 29, 57, 84)


