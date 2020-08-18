library(ape)
misof <- read.nexus("misof.backbone.nex")
max(branching.times(misof[[1]]))
sum(misof[[1]]$edge.length)

rain <- read.nexus("rainford.backbone.nex")
max(branching.times(rain[[1]]))
sum(rain[[1]]$edge.length)


sum(rain[[1]]$edge.length)/sum(misof[[1]]$edge.length)
