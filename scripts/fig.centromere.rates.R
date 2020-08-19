# this script makes the plot comparing delta R stats for
# monocentric and holocentric
par(mfcol=c(1,2))
load("../results/cent.rates-backbone1.RData")
library(coda)
fission <- rates$asc1-rates$asc2
fusion <- rates$desc1-rates$desc2
poly <- rates$pol1-rates$pol2
hpdfis <- HPDinterval(as.mcmc(fission))
hpdfus <- HPDinterval(as.mcmc(fusion))
hpdpol <- HPDinterval(as.mcmc(poly))

#state 1 = holocentric, state 2 = monocentric
HPDinterval(as.mcmc(rates$asc1))
HPDinterval(as.mcmc(rates$asc2))

HPDinterval(as.mcmc(rates$desc1))
HPDinterval(as.mcmc(rates$desc2))

HPDinterval(as.mcmc(rates$pol1))
HPDinterval(as.mcmc(rates$pol2))

colMeans(rates)


cols <- c(rgb(1, 0, 0, .5), rgb(0, 1, 0, .5), rgb(0, 0, 1, .5))
plot(0,0,col="white",
     ylim=c(-50,700),
     xlim=c(-.02,.01),
     xlab=expression(paste(Delta, R[x])),
     ylab="density")
abline(v=0, lty=3,col="black")
polygon(density(fission, bw = .0005), col = cols[1])
polygon(density(fusion, bw = .0005),col = cols[2])
polygon(density(poly, bw = .0005),col = cols[3])
points(pch = 22, bg = cols,
       x = rep(-.02, 3), y = c(700, 650, 600))
text(x = rep(-.02, 3), y = c(700, 650, 600), pos = 4, cex = .7,
     labels=c("fission", "fusion", "polyploidy"))

lines(y=rep(-15, 2), x=hpdfis, col=cols[1], lwd=4)
lines(y=rep(-37, 2), x=hpdfus, col=cols[2], lwd=4)
lines(y=rep(-59, 2), x=hpdpol, col=cols[3], lwd=4)
# export this as pdf at 8x4



# now for the simple model

load("../results/cent.rates.simple.model.RData")
#rates must be back transformed into units of
# millions of years tree depth is 479.1039
rates[,2:7] <- rates[,2:7]/479.1039

library(coda)
fission <- rates$asc1-rates$asc2
fusion <- rates$desc1-rates$desc2
hpdfis <- HPDinterval(as.mcmc(fission))
hpdfus <- HPDinterval(as.mcmc(fusion))


#state 1 = holocentric, state 2 = monocentric
HPDinterval(as.mcmc(rates$asc1))
HPDinterval(as.mcmc(rates$asc2))

HPDinterval(as.mcmc(rates$desc1))
HPDinterval(as.mcmc(rates$desc2))

colMeans(rates)


cols <- c(rgb(1, 0, 0, .5), rgb(0, 1, 0, .5))
plot(0,0,col="white",
     ylim=c(-5,75),
     xlim=c(-.06,.1),
     xlab=expression(paste(Delta, R[x])),
     ylab="density")
abline(v=0, lty=3,col="black")
polygon(density(fission), col = cols[1])
polygon(density(fusion),col = cols[2])
lines(y=rep(-1.7, 2), x=hpdfis, col=cols[1], lwd=4)
lines(y=rep(-4, 2), x=hpdfus, col=cols[2], lwd=4)
# export this as pdf at 4x4

