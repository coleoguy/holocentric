# this script contains chunks of code to create our
# plot of our data these pieces are then combined for
# the final figure presented in the paper.

library(phytools)
library(chromePlus)
library(viridis)
trees <- read.nexus("../data/misof.backbone.nex")
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

# export 5x5

# get bar tips
chrom.num <- dat$haploid
names(chrom.num) <- dat$genus
foo <- read.csv("../data/phylo.data.csv", as.is=T)
tip.orders <- c()
for(i in 1:599){
  tip.orders[i] <- foo$Order[foo$Genus == tree$tip.label[i]][1]
}
tip.orders[tip.orders=="Hempitera"] <- "Hemiptera"
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
               border=NA,
               width=.009,
               edge.width=.01)
# export at 5x5


