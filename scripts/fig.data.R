library(phytools)
library(chromePlus)
library(viridis)
trees <- read.nexus("../data/misof.backbone.nex")
#trees <- read.nexus("../data/rainford.backbone.nex")
dat <- read.csv("../data/data.invert.csv", as.is=T)


source("functions.R")

foo <- getData2(trees, dat)
tree <- foo[[1]][[1]]
dat <- foo[[2]]
# mono is zero

chrom.type <- dat$chrom
names(chrom.type) <- dat$genus
new.haps <- round(log(dat$haploid)*10)
range(new.haps)



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










# get bar tips
chrom.num <- dat$haploid
names(chrom.num) <- dat$genus
foo <- read.csv("../data/phylo.data.csv", as.is=T)
tip.orders <- c()
for(i in 1:599){
  tip.orders[i] <- foo$Order[foo$Genus == tree$tip.label[i]][1]
}
tip.colors <- rep("darkgray", 599)
tip.colors[tip.orders == "Lepidoptera"] <- "#fa7763"
tip.colors[tip.orders == "Hemiptera"] <- "#e59404"
tip.colors[tip.orders == "Odonata"] <- "#a2a903"
tip.colors[tip.orders == "Blattodea"] <- "#28b600"
tip.colors[tip.orders == "Phasmatodea"] <- "#00c07b"
tip.colors[tip.orders == "Isoptera"] <- "#00c4c3"
tip.colors[tip.orders == "Coleoptera"] <- "#0ab6ff"
tip.colors[tip.orders == "Hymenoptera"] <- "#9892ff"
tip.colors[tip.orders == "Neuroptera"] <- "#e56cf5"
tip.colors[tip.orders == "Diptera"] <- "#ff5fc0"

plotTree.wBars(tree=tree,
               x=chrom.num,
               type="fan",
               col=tip.colors,
               border=NA)



