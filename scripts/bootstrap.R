#read in data
rates <- load("../results/cent.rates-backbone1.RData")


starts <- seq(from=1, by=25, length.out = 100)
ends <- seq(from=25, by=25, length.out = 100)
treerows <- matrix(, 25, 100)
for(i in 1:100){
  treerows[,i] <- starts[i]:ends[i]
}
sig.count <- 0

library(coda)
dists <- list()
for(i in 1:1000){
  sel.trees <- sample(1:100, replace=T)
  sel.rows <- as.vector(treerows[,sel.trees])
  if(HPDinterval(as.mcmc(rates$asc1[sel.rows]))[2] <
     HPDinterval(as.mcmc(rates$asc2[sel.rows]))[1]){
    sig.count <- sig.count + 1
  }
 dists[[i]] <- matrix(c(rates$asc1[sel.rows], rates$asc2[sel.rows]),
                      2500, 2)
}

plot(0,0,col="white",
     ylim=c(-5,110),
     xlim=c(-.05,.1),
     xlab=expression(paste(Delta, R[x])),
     ylab="density")
abline(v=0, lty=3,col="black")
cols <- rainbow(1000)
for(i in 1:1000){
rtemp <- dists[[i]][,1] + dists[[i]][,2]

lines(density(rtemp, bw = .0005), col = cols[1])}
polygon(density(fusion, bw = .0005),col = cols[2])
polygon(density(poly, bw = .0005),col = cols[3])
points(pch = 22, bg = cols,
       x = rep(-.02, 3), y = c(700, 650, 600))
text(x = rep(-.02, 3), y = c(700, 650, 600), pos = 4, cex = .7,
     labels=c("fission", "fusion", "polyploidy"))

lines(y=rep(-15, 2), x=hpdfis, col=cols[1], lwd=4)
lines(y=rep(-37, 2), x=hpdfus, col=cols[2], lwd=4)
lines(y=rep(-59, 2), x=hpdpol, col=cols[3], lwd=4)






abline(v=starts, col="red", lwd=.3)
lines(rates$asc2, col="blue")
library(coda)
HPDinterval(as.mcmc(rates$asc1[1:25]))
HPDinterval(as.mcmc(rates$asc2[1:25]))

plot(rates$asc1~rates$asc2, pch=16, cex=.3)
summary(lm(rates$asc1~rates$asc2))
