load("../results/cent.rates.RData")

plot(0,0,col="white",
     ylim=c(0,400),
     xlim=c(-.02,.02),
     xlab="rate difference",
     ylab="density")
abline(v=0, lty=3,col="red")
polygon(density(rates$asc1-rates$asc2,bw=.001), col=rgb(1,0,0,.1))
polygon(density(rates$desc1-rates$desc2,bw=.001),col=rgb(0,1,0,.1))
polygon(density(rates$pol1-rates$pol2,bw=.001),col=rgb(0,0,1,.1))
points(pch=22, bg=c(rgb(1,0,0,.2), rgb(0,1,0,.2), rgb(0,0,1,.2)),
       x=rep(.01,3), y=c(400,370,340))
text(x=rep(.01,3), y=c(400,370,340), pos=4, cex=.7,labels=c("ascending","descending","polyploidy"))
# export this as pdf at 4x4
