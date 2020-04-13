# this function takes posterior of trees and data
# it returns a list with pruned data and random sample of possible
# tip values if multiple points are available.
getData <- function(trees, dat){
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
  tree.depths <- c()
  for(i in 1:100){
    cur.tree <- drop.tip(trees[[i]], tip = missing)
    cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
    tree.depths[i] <- max(branching.times(cur.tree))
    trees.pruned[[i]] <- cur.tree
  }
  chrom <- data.frame(dat.pruned$Genus,
                      as.numeric(dat.pruned$haploid.num),
                      dat.pruned$chromosome, stringsAsFactors = F)
  colnames(chrom) <- c("genus", "haploid", "chrom")
  # this code means that monocentric will be state 2 in the output
  # and holocentric will be state one in the output
  chrom$chrom[chrom$chrom == "mono"] <- 0
  chrom$chrom[chrom$chrom == "holo"] <- 1
  #chrom <- chrom[complete.cases(chrom), ]
  chrom$chrom <- as.numeric(chrom$chrom)
  chroms <- datatoMatrix(chrom,
                         range = range(chrom$haploid) + c(0, 3),
                         hyper = T)

  results <- list()
  results[[1]] <- trees.pruned
  results[[2]] <- chroms
  results[[3]] <- tree.depths
  return(results)
}












# this function takes posterior of trees and data
# it returns a list with pruned data and random sample of possible
# tip values if multiple points are available.
getData2 <- function(trees, dat){
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
  tree.depths <- c()
  for(i in 1:100){
    cur.tree <- drop.tip(trees[[i]], tip = missing)
    cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
    tree.depths[i] <- max(branching.times(cur.tree))
    trees.pruned[[i]] <- cur.tree
  }
  chrom <- data.frame(dat.pruned$Genus,
                      as.numeric(dat.pruned$haploid.num),
                      dat.pruned$chromosome, stringsAsFactors = F)
  colnames(chrom) <- c("genus", "haploid", "chrom")
  # this code means that monocentric will be state 2 in the output
  # and holocentric will be state one in the output
#  chrom$chrom[chrom$chrom == "mono"] <- 0
#  chrom$chrom[chrom$chrom == "holo"] <- 1
  #chrom <- chrom[complete.cases(chrom), ]
  #chrom$chrom <- as.numeric(chrom$chrom)
  #chroms <- datatoMatrix(chrom,
  #                       range = range(chrom$haploid) + c(0, 3),
   #                      hyper = T)

  results <- list()
  results[[1]] <- trees.pruned
  results[[2]] <- chrom
  results[[3]] <- tree.depths
  return(results)
}


# this function takes posterior of trees and data
# it returns a list with pruned data and random sample of possible
# tip values if multiple points are available.
getDataOrder <- function(trees, dat, order){
  tree.genera <- trees[[1]]$tip.label
  good.genera <- unique(dat$Genus[which(dat$Genus %in% tree.genera)])
  good.genera <- good.genera[good.genera %in% dat$Genus[dat$Order==order]]
  hit <- sample(which(dat$Genus == good.genera[1]), 1)
  dat.pruned <- dat[hit, ]
  for(i in 2:length(good.genera)){
    hit <- which(dat$Genus == good.genera[i])
    if(length(hit)>1)  hit <- sample(hit, 1)
    dat.pruned <- rbind(dat.pruned, dat[hit, ])
  }
  missing <- trees[[1]]$tip.label[!trees[[1]]$tip.label %in% dat.pruned$Genus]
  trees.pruned <- list()
  tree.depths <- c()
  for(i in 1:100){
    cur.tree <- drop.tip(trees[[i]], tip = missing)
    cur.tree$edge.length <-  cur.tree$edge.length / max(branching.times(cur.tree))
    tree.depths[i] <- max(branching.times(cur.tree))
    trees.pruned[[i]] <- cur.tree
  }
  chrom <- data.frame(dat.pruned$Genus,
                      as.numeric(dat.pruned$haploid.num),
                      dat.pruned$chromosome, stringsAsFactors = F)
  colnames(chrom) <- c("genus", "haploid", "chrom")
  #this code means that monocentric will be state 2 in the output
  #and holocentric will be state one in the output
  chrom$chrom[chrom$chrom == "mono"] <- 0
  chrom$chrom[chrom$chrom == "holo"] <- 1
  chrom <- chrom[complete.cases(chrom), ]
  chrom$chrom <- as.numeric(chrom$chrom)
  chroms <- datatoMatrix(chrom,
                         range = range(chrom$haploid) + c(0, 3),
                        hyper = F)

  results <- list()
  results[[1]] <- trees.pruned
  results[[2]] <- chroms
  results[[3]] <- tree.depths
  return(results)
}

