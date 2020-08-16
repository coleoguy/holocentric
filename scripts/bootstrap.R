#read in data
rates <- load("../results/cent.rates-backbone1.RData")

starts <- seq(from=1, by=25, length.out = 100)
ends <- seq(from=25, by=25, length.out = 100)
treerows <- matrix(, 25, 100)
for(i in 1:100){
  treerows[,i] <- starts[i]:ends[i]
}


# MAKE SUP PLOTS

# get fissions
library(coda)
dists <- matrix(NA,2500,1000)
for(i in 1:1000){
  sel.trees <- sample(1:100, replace=T)
  sel.rows <- as.vector(treerows[,sel.trees])
  dists[, i] <- rates$asc1[sel.rows] - rates$asc2[sel.rows]
}
par(mfcol=c(1,3))
# plot fissions
library(viridis)
plot(0,0,col="white",
     ylim=c(-5,105),
     xlim=c(-.025,.04),
     xlab=expression(paste(Delta, R[fission])),
     ylab="density")
cols <- sample(viridis(1000))
lows <- highs <- c()
for(i in 1:1000){
  lines(density(dists[, i], bw = .0015), col = cols[i],lwd=.09)
  highs[i] <- HPDinterval(as.mcmc(dists[,i]))[2]
  lows[i] <- HPDinterval(as.mcmc(dists[,i]))[1]
}
lines(density(rates$asc1 - rates$asc2,, bw = .0015), lty=3, lwd=2)
lines(y=c(-5,-5),x=c(max(lows),min(highs)),lwd=3)

# get fusions
dists <- matrix(NA,2500,1000)
for(i in 1:1000){
  sel.trees <- sample(1:100, replace=T)
  sel.rows <- as.vector(treerows[,sel.trees])
  dists[, i] <- rates$desc1[sel.rows] - rates$desc2[sel.rows]
}

# plot fusions
plot(0,0,col="white",
     ylim=c(-5,165),
     xlim=c(-.025,.03),
     xlab=expression(paste(Delta, R[fusion])),
     ylab="density")
for(i in 1:1000){
  lines(density(dists[, i], bw = .0015), col = cols[i],lwd=.09)
  highs[i] <- HPDinterval(as.mcmc(dists[,i]))[2]
  lows[i] <- HPDinterval(as.mcmc(dists[,i]))[1]
}
lines(density(rates$desc1 - rates$desc2,, bw = .0015), lty=3, lwd=2)
lines(y=c(-5,-5),x=c(max(lows),min(highs)),lwd=3)

# get polyploidy
dists <- matrix(NA,2500,1000)
for(i in 1:1000){
  sel.trees <- sample(1:100, replace=T)
  sel.rows <- as.vector(treerows[,sel.trees])
  dists[, i] <- rates$pol1[sel.rows] - rates$pol2[sel.rows]
}

# plot polyploidy
plot(0,0,col="white",
     ylim=c(-5,1400),
     xlim=c(-.002,.002),
     xlab=expression(paste(Delta, R[polyploidy])),
     ylab="density")
for(i in 1:1000){
  lines(density(dists[, i], bw = .0001), col = cols[i],lwd=.09)
  highs[i] <- HPDinterval(as.mcmc(dists[,i]))[2]
  lows[i] <- HPDinterval(as.mcmc(dists[,i]))[1]
}
lines(density(rates$pol1 - rates$pol2,, bw = .0001), lty=3, lwd=2)
lines(y=c(-5,-5),x=c(max(lows),min(highs)),lwd=3)

# export 10" wide by 4" tall


