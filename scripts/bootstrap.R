rates <- f.results
starts <- seq(from=1, by=25, length.out = 100)
ends <- seq(from=25, by=25, length.out = 100)
treerows <- matrix(, 25, 100)
for(i in 1:100){
  treerows[,i] <- starts[i]:ends[i]
}
sig.count <- 0
for(i in 1:1000){
  sel.trees <- sample(1:100, replace=T)
  sel.rows <- as.vector(treerows[,sel.trees])
  if(HPDinterval(as.mcmc(rates$asc1[sel.rows]))[2] <
     HPDinterval(as.mcmc(rates$asc2[sel.rows]))[1]){
    sig.count <- sig.count + 1
  }
}







abline(v=starts, col="red", lwd=.3)
lines(rates$asc2, col="blue")
library(coda)
HPDinterval(as.mcmc(rates$asc1[1:25]))
HPDinterval(as.mcmc(rates$asc2[1:25]))

plot(rates$asc1~rates$asc2, pch=16, cex=.3)
summary(lm(rates$asc1~rates$asc2))
