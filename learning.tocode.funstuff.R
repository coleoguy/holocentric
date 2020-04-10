# Heath, Michelle, and Sarah testing this old holocentric
# monocentric idea on a tree of insects

# load packages
library(ape)
library(chromePlus)
library(diversitree)

# load custom functions
source("functions.R")

# read in the data and the tree
trees <- read.nexus("misof.backbone.nex")
dat <- read.csv("data.invert.csv", as.is = T)[, -c(6:9)]

# set the MCMC chain length
iter <- 10

# we scale our trees so lets store that info for transforming
# rates back into millions of years
rate.scalers <- c()

# we will loop through all 100 trees
for(i in 1:2){
  # slim the data to include only the desired data
  # and generate the format table needed by chromPlus
  foo <- getData(trees, dat)
  trees.pruned <- foo[[1]]
  chroms <- foo[[2]]
  rate.scalers[i] <- foo[[3]][i]
  rm(foo)

  # make the basic likelihood function for the data
  lk.mk <- make.mkn(trees.pruned[[i]], states = chroms,
                    k = ncol(chroms), strict = F,
                    control = list(method = "ode"))

  # now we constrain our model to be biologically realistic for
  # chromosomes.
  con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = T,
                          polyploidy = F, verbose = F,
                          constrain = list(drop.demi = T, drop.poly = T))

  # this loop does a test run to identify good w values
  # and sets up priors we only need to run it one time
  if(i == 1){
    argnames(con.lk.mk)
    prior <- make.prior.exponential(1)
    temp <- mcmc(con.lk.mk,
                 x.init = runif(6, 0, 10),
                 prior = prior,
                 w = 1,
                 nsteps = iter/10)
    temp <- temp[-c(1:10), ]
    w <- diff(sapply(temp[2:7],
                     quantile, c(.05, .95)))
  }

  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init = as.numeric(temp[nrow(temp), 2:7]),
                      prior = prior,
                      w = w, nsteps = iter)

  # just in case we have a crash lets write results for each tree
  write.csv(results[[i]], file=paste("tree",i,".csv", sep=""))
}

# NOTE monocentric is state 2 and holocentric is state 1

