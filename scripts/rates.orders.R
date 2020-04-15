# Heath, Michelle, and Sarah testing this old holocentric
# monocentric idea on a tree of insects

# load packages
library(ape) # basic phylo tools
library(chromePlus) # chromosome models
library(diversitree) # basic likelihood functions
library(doMC) # this allows multicore runs on a mac
registerDoMC(14)
iter <- 50
# for troubleshooting we might run on just a couple of trees
# usually ntree should be set equal to the number of trees being
# analyzed.
ntree <- 100

# load custom functions
source("functions.R")

# read in the data and the tree
trees <- read.nexus("../data/misof.backbone.nex")
dat <- read.csv("../data/data.invert.csv", as.is = T)[, -c(6:9)]

# lets find our orders first
foo <- getData(trees, dat)
tree <- foo[[1]]
genera <- row.names(foo[[2]])
orders <- dat$Order[dat$Genus %in% genera]
orders <- as.data.frame(table(orders))
# here I am reducing down to only estimate rates in orders
# with at least 20 species
orders <- as.character(orders$orders[orders$Freq>=20])

# this large outside loop will repeat for each order we are
# evaluating.

for(i in 8:length(orders)){
  foo <- getDataOrder(trees, dat, order = orders[i])
  trees.pruned <- foo[[1]]
  chroms <- foo[[2]]
  lk.mk <- make.mkn(trees.pruned[[1]], states = chroms,
                    k = ncol(chroms), strict = F,
                    control = list(method = "ode"))
  # we will do this next bit twice to get w with and without polyploidy
  # first with polyploidy
  con.lk.mk.wop<-constrainMkn(data = chroms, lik = lk.mk, hyper = F,
                              polyploidy = F, verbose = F,
                              constrain = list(drop.demi = T, drop.poly = T))
  con.lk.mk.wp<-constrainMkn(data = chroms, lik = lk.mk, hyper = F,
                             polyploidy = F, verbose = F,
                             constrain = list(drop.demi = T, drop.poly = F))
  prior <- make.prior.exponential(.5)
  temp.wop <- mcmc(con.lk.mk.wop,
                   x.init = runif(2, 0, 1),
                   prior = prior,
                   w = 1,
                   nsteps = 40,
                   upper = 50,
                   lower = 0)
  temp.wp <- mcmc(con.lk.mk.wp,
                   x.init = runif(3, 0, 1),
                   prior = prior,
                   w = 1,
                   nsteps = 40,
                   upper = 50,
                   lower = 0)
  temp.wop <- temp.wop[-c(1:10), ]
  temp.wp <- temp.wp[-c(1:10), ]
  w.wop <- diff(sapply(temp.wop[2:3],
                      quantile, c(.05, .95)))
  w.wp <- diff(sapply(temp.wp[2:4],
                      quantile, c(.05, .95)))
  x <- foreach(j=1:ntree) %dopar%{
    # slim the data to include only the desired data
    # and generate the format table needed by chromPlus
    foo <- getDataOrder(trees, dat, order=orders[i])
    trees.pruned <- foo[[1]]
    chroms <- foo[[2]]
    rm(foo)
    lk.mk <- make.mkn(trees.pruned[[j]], states = chroms,
                      k = ncol(chroms), strict = F,
                      control = list(method = "ode"))
    con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = F,
                            polyploidy = F, verbose = F,
                            constrain = list(drop.demi = T, drop.poly = T))
    cur.res <- mcmc(con.lk.mk,
                        x.init =  runif(2, 0, 1),
                        prior = prior,
                        w = w.wop,
                        nsteps = iter,
                        upper = 50,
                        lower = 0)
    cur.res
  }
  scaler <- getDataOrder(trees, dat, order=orders[i])[[3]]
  for(k in 1:ntree){
    x[[k]][,2:3] <- x[[k]][,2:3]/scaler[k]
  }
  save(x, file=paste("../results/rates.", orders[i], ".wop.rda", sep=""))

#### NOW WE RUN WITH POLYPLOIDY
  x <- foreach(j=1:ntree) %dopar%{
    # slim the data to include only the desired data
    # and generate the format table needed by chromPlus
    foo <- getDataOrder(trees, dat, order=orders[i])
    trees.pruned <- foo[[1]]
    chroms <- foo[[2]]
    rm(foo)
    lk.mk <- make.mkn(trees.pruned[[j]], states = chroms,
                      k = ncol(chroms), strict = F,
                      control = list(method = "ode"))
    con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = F,
                            polyploidy = F, verbose = F,
                            constrain = list(drop.demi = T, drop.poly = F))
    cur.res <- mcmc(con.lk.mk,
                    x.init =  runif(3, 0, 1),
                    prior = prior,
                    w = w.wp,
                    nsteps = iter,
                    upper = 50,
                    lower = 0)
    cur.res
  }
  scaler <- getDataOrder(trees, dat, order=orders[i])[[3]]
  for(k in 1:ntree){
    x[[k]][,2:4] <- x[[k]][,2:4]/scaler[k]
  }
  save(x, file=paste("../results/rates.", orders[i], ".wp.rda", sep=""))
}





















