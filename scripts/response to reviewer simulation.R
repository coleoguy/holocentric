# load packages
library(chromePlus)
library(diversitree)

# set seed
set.seed(1)

# simulate a phylogeny under birth death model
tree <- tree.bd(pars = c(3,1), max.taxa=100)

# scale the tree to unit length
tree$edge.length <- tree$edge.length/max(branching.times(tree))

# simulate chromosome number evolution
# we will simulate 100 datasets with
# a fusion and fission rate of .1 (slowdat)
# and an additional 100 datasets with
# a fusion and fission rate of 1 (fastdat)
# for each simulation we will record the minimum and maximum
# chromosome number
fastmax <- fastmin <- slowmax <- slowmin <- c()
for(i in 1:1000){
  print(i)
  temp <- simChrom(tree, pars=c(1, 1, 0, 0, 10),
                           limits = c(1, 20), model = "2010")
  fastmax[i] <- max(temp)
  fastmin[i] <- min(temp)
  temp <- simChrom(tree, pars=c(.1, .1, 0, 0, 10),
                           limits = c(1, 20), model = "2010")
  slowmax[i] <- max(temp)
  slowmin[i] <- min(temp)
}
df <- data.frame(c(fastmax, fastmin, slowmax, slowmin),
                 rep(c("fastmax","fastmin","slowmax","slowmin"),
                     each=1000))
colnames(df) <- c("chromnumber","model")
ggplot(data = df, aes(x = chromnumber, fill = model)) +
  geom_histogram(position = "dodge", binwidth = 1)
