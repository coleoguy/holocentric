# Heath, Michelle, and Sarah testing this old holocentric
# monocentric idea on a tree of insects

# load packages
library(ape) # basic phylo tools
library(chromePlus) # chromosome models
library(diversitree) # basic likelihood functions
library(doMC) # this allows multicore runs on a mac
registerDoMC(14)

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
orders <- as.character(orders$orders[orders$Freq>=20])

# this large outside loop will repeat for each order we are
# evaluating.

for(i in 1:length(orders)){
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
                   nsteps = 20,
                   upper = 50,
                   lower = 0)
  temp.wp <- mcmc(con.lk.mk.wp,
                   x.init = runif(3, 0, 1),
                   prior = prior,
                   w = 1,
                   nsteps = 20,
                   upper = 50,
                   lower = 0)
  temp.wop <- temp.wop[-c(1:10), ]
  temp.wp <- temp.wp[-c(1:10), ]
  w.wop <- diff(sapply(temp.wop[2:3],
                      quantile, c(.05, .95)))
  w.wp <- diff(sapply(temp.wp[2:4],
                      quantile, c(.05, .95)))
  ## OK now we are ready to run our analysis
  results <- list()
  x <- foreach(j=1:100) %dopar%{
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
    result[[j]] <- mcmc(con.lk.mk,
                        x.init =  runif(2, 0, 1),
                        prior = prior,
                        w = w.wop,
                        nsteps = iter,
                        upper = 50,
                        lower = 0)
    write.csv(result[[j]], file=paste("tree.nop",
                                      orders[i],
                                      ".", j,".csv",
                                      sep=""))
  }
  results <- list()
  for(i in 1:100){
    scaler <- max(branching.times(trees.pruned[[i]]))
    x[,2:3] <- x[,2:3]/scaler
    results[[i]] <- x
  }
  save(results, file="rates.", orders[i] ,".wop")

  results <- list()
  x <- foreach(j=1:100) %dopar%{
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
    result[[j]] <- mcmc(con.lk.mk,
                        x.init =  runif(3, 0, 1),
                        prior = prior,
                        w = w.wp,
                        nsteps = iter,
                        upper = 50,
                        lower = 0)
    write.csv(result[[j]], file=paste("tree.p",
                                      orders[i],
                                      ".", j,".csv",
                                      sep=""))
  }
  results <- list()
  for(i in 1:100){
    scaler <- max(branching.times(trees.pruned[[i]]))
    x[,2:4] <- x[,2:4]/scaler
    results[[i]] <- x
  }
  save(results, file="rates.", orders[i] ,".wp")






}















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
rm(foo)



result <- list()
# we will loop through all 100 trees

# NOTE monocentric is state 2 and holocentric is state 1 in results

for(i in 1:10){
  foo <- read.csv(paste("tree", i,".csv", sep=""))
  if(i == 1){
    plot(foo$pol1 - foo$pol2,type="l", ylim=c(-.4,.6))
  }else{
    lines(foo$pol1 - foo$pol2, col=rainbow(10)[i])
  }
}



