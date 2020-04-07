# Heath pretending to be Michelle or Sarah

library(ape)
# read in the data and the tree
trees <- read.nexus("misof.backbone.nex")
dat <- read.csv("data.invert.csv", as.is=T)[,-c(6:9)]


# TODO compare tree and data and get just the overlap
tree.genera <- trees[[1]]$tip.label
good.genera <- unique(dat$Genus[which(dat$Genus %in% tree.genera)])
hit <- sample(which(dat$Genus == good.genera[1]), 1)
dat.pruned <- dat[hit, ]
for(i in 2:length(good.genera)){
  hit <- which(dat$Genus == good.genera[i])
  if(length(hit)>1)  hit <- sample(hit, 1)
  dat.pruned <- rbind(dat.pruned, dat[hit, ])
}

missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% dat.pruned$Genus]
trees.pruned <- list()
for(i in 1:100){
  cur.tree <- drop.tip(trees[[i]], tip = missing)
  cur.tree$edge.length <-  cur.tree$edge.length/ max(branching.times(cur.tree))
  trees.pruned[[i]] <- cur.tree
}


library(chromePlus)

# slim the data to include only the desired data
chrom <- data.frame(dat.pruned$Genus,
                    dat.pruned$haploid.num,
                    dat.pruned$chromosome, stringsAsFactors = F)
colnames(chrom) <- c("genus", "haploid", "chrom")
chrom$chrom[chrom$chrom == "mono"] <- 0
chrom$chrom[chrom$chrom == "holo"] <- 1
chrom <- chrom[complete.cases(chrom), ]
chrom$haploid <- as.numeric(chrom$haploid)
chrom$chrom <- as.numeric(chrom$chrom)
### range should be larger than that observed in the
### data however it should be as small as possible without
### impacting inference to speed computation
# check the range of chromosome numbers present
range(chrom$haploid)

chroms <- datatoMatrix(chrom, range=range(chrom$haploid), hyper=T)

# set the MCMC chain length
iter <- 20
library(diversitree)

### w is best set by doing a short run and getting some idea of the
### confidence interval on these parameters
lk.mk <- make.mkn(trees.pruned[[1]], states=chroms, k=ncol(chroms), strict=F, control=list(method="ode"))

### This sets up the model described in the paper
con.lk.mk<-constrainMkn(data=chroms, lik=lk.mk, hyper=T,
                        polyploidy=F, verbose=F,
                        constrain=list(drop.demi=T))
argnames(con.lk.mk)
temp <- mcmc(con.lk.mk, x.init=runif(8,0,10), w=1, nsteps=iter)
temp[-c(1:5),]
w <- diff(sapply(temp[2:9], quantile, c(.05, .95)))


### this will cycle through all 100 trees
results <- list()
for(i in 1:100){
  # make the initial likelihood function
  lk.mk <- make.mkn(trees.pruned[[i]], states = chroms, k = ncol(chroms),
                    strict = F, control = list(method = "ode"))
  # constrain to a biologically realistic model of chrom evolution
  con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T))
  # run the MCMC
  results[[i]] <- mcmc(con.lk.mk, x.init = colMeans(temp)[2:9],
                       w = w, nsteps = iter)
}


