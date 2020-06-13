# This script make the supplemental figure that shows the impact
# of using different backbone trees

load("../results/cent.rates-backbone2.RData")
load("../results/cent.rates-backbone1.RData")
ratesb1 <- rates[,-1]
ratesb2 <- f.results
rm(f.results, rates)


par(mfcol=c(2,2))
plot(density(ratesb1$asc1 - ratesb1$asc2), ylim=c(0,120),
     xlab=expression(paste(Delta, R[gamma])),
     main="Fission")
lines(density(ratesb2$asc1 - ratesb2$asc2),col="red",
      main="Fission")
abline(v=0,lty=2)

plot(density(ratesb1$desc1 - ratesb1$desc2), ylim=c(0,240),
     xlab=expression(paste(Delta, R[delta])),
     main="Fusion")
lines(density(ratesb2$desc1 - ratesb2$desc2),col="red")
abline(v=0,lty=2)

plot(density(ratesb1$pol1 - ratesb1$pol2), ylim=c(0,1520),
     main="Polyploidy",
     xlab=expression(paste(Delta, R[rho])))
lines(density(ratesb2$pol1 - ratesb2$pol2),col="red")
abline(v=0,lty=2)

