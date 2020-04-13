# Heath, Michelle, and Sarah testing this old holocentric
# monocentric idea on a tree of insects

# load packages
library(ape) # basic phylo tools
library(chromePlus) # chromosome models
library(diversitree) # basic likelihood functions
library(doMC) # this allows multicore runs on a mac

# load custom functions
source("functions.R")

# read in the data and the tree
trees <- read.nexus("misof.backbone.nex")
dat <- read.csv("data.invert.csv", as.is = T)[, -c(6:9)]

# set the MCMC chain length
iter <- 50

# we scale our trees so lets store that info for transforming
# rates back into millions of years
rate.scalers <- c()

#############################
#                           #
# TEST RUN TO GET W         #
#                           #
#############################

# slim the data to include only the desired data
# and generate the format table needed by chromPlus
foo <- getData(trees, dat)
trees.pruned <- foo[[1]]
chroms <- foo[[2]]
rm(foo)
lk.mk <- make.mkn(trees.pruned[[1]], states = chroms,
                  k = ncol(chroms), strict = F,
                  control = list(method = "ode"))
con.lk.mk<-constrainMkn(data = chroms, lik = lk.mk, hyper = T,
                        polyploidy = F, verbose = F,
                        constrain = list(drop.demi = T, drop.poly = T))
prior <- make.prior.exponential(.5)
temp <- mcmc(con.lk.mk,
             x.init = runif(6, 0, 1),
             prior = prior,
             w = 1,
             nsteps = 20,
             upper = 25,
             lower = 0)
temp <- temp[-c(1:10), ]
w <- diff(sapply(temp[2:7],
                 quantile, c(.05, .95)))

#################################
#                               #
#  Now we can do our full run   #
#                               #
#################################

# this will allow us to run on 12 cores
registerDoMC(10)
result <- list()
# we will loop through all 100 trees
x <- foreach(i=1:10) %dopar%{
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

  # now we are ready to run our inference run
  result[[i]] <- mcmc(con.lk.mk,
                      x.init =  runif(6, 0, 1),
                      prior = prior,
                      w = w,
                      nsteps = iter,
                      upper = 25,
                      lower = 0)

  # just in case we have a crash lets write results for each tree
  write.csv(result[[i]], file=paste("tree.nop",i,".csv", sep=""))
}

# NOTE monocentric is state 2 and holocentric is state 1 in results

for(i in 1:10){
  foo <- read.csv(paste("tree", i,".csv", sep=""))
  if(i == 1){
    plot(foo$pol1 - foo$pol2,type="l", ylim=c(-.4,.6))
  }else{
    lines(foo$pol1 - foo$pol2, col=rainbow(10)[i])
  }
}



