library(ggplot2)
library(tidybayes)
load(file="../results/order.rates.RData")

# get parameter names correct
foo <- as.character(results.wop$type)
foo[foo=="asc"] <- "fission"
foo[foo=="desc"] <- "fusion"
results.wop$type <- as.factor(foo)

foo <- as.character(results.wp$type)
foo[foo=="asc"] <- "fission"
foo[foo=="pol"] <- "polyploidy"
foo[foo=="desc"] <- "fusion"
results.wp$type <- as.factor(foo)

# lets figure out the order for the data
foo <- aggregate(x = results.wp$rate,
                 by = list(results.wp$order, results.wp$type),
                 FUN = mean)
foo <- foo[foo$Group.2=="fission",]
foo <- foo[order(foo$x, decreasing=T),]
x <- row.names(foo)[which(foo$Group.1 %in% c("Hemiptera","Lepidoptera","Odonata"))]
x <- c(x, row.names(foo)[which(!foo$Group.1 %in% c("Hemiptera","Lepidoptera","Odonata"))])
x <- as.numeric(x)

# reorder so taxa are in right order
results.wop$order <- factor(results.wop$order,
                            levels(results.wop$order)[x])
results.wp$order <- factor(results.wp$order,
                           levels(results.wp$order)[x])


# this code is for single legend
ggplot(results.wop, aes(x=type, y=rate, color=order)) +
  geom_jitter(cex=2, alpha=1,position=position_jitterdodge(.3)) +
  theme_bw()+
  guides(color=guide_legend(title="Order"))
# this code is for plot wo polyploidy
ggplot(results.wop, aes(x = type, y = rate, color = order)) +
  geom_jitter(cex = .5, alpha = .1,
              position = position_jitterdodge(.3)) +
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Parameter") + ylab("Rate") +
  stat_summary(aes(x = type, y = rate, fill = order),
               fun.data = "mean_hdci", fun.args = list(mult=1),
               size = 0.4, position = position_jitterdodge(0),
               inherit.aes = FALSE) # exported 4x4
# this code is for plot w polyploidy
ggplot(results.wp, aes(x=type, y=rate, color=order)) +
  geom_jitter(cex=.5, alpha=.1,position=position_jitterdodge(.3)) +
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Parameter") + ylab("Rate") +
  stat_summary(aes(x=type, y=rate, fill=order), fun.data="mean_hdci", fun.args = list(mult=1),
               size = 0.4, position = position_jitterdodge(0),
               inherit.aes = FALSE) # exported 6x4

