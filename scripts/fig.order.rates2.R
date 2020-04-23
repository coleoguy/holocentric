library(ggraptR)
library(ggplot2)
#ggraptR(results.wop)


load(file="../results/order.rates.RData")

results.wop$order <- factor(results.wop$order, levels(results.wop$order)[c(4, 7, 9, 1:3, 5, 6, 8, 10)])
results.wp$order <- factor(results.wp$order, levels(results.wp$order)[c(4, 7, 9, 1:3, 5, 6, 8, 10)])

#Plot of results without polyploidy
ggplot(results.wop, aes(y=rate, x=as.factor(type))) +
  geom_boxplot(aes(fill=as.factor(order)), stat="boxplot",
               alpha=0.5, width=0.4, position=position_dodge(0.7)) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000",
                          size=15, hjust=0.5, vjust=0.5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  guides(fill=guide_legend(title="Order")) +
   xlab("Type of Rearrangement") + ylab("Rate")
# export at 7"x3.5"
#Plot of results with polyploidy
ggplot(results.wp, aes(y=rate, x=as.factor(type))) +
  geom_boxplot(aes(fill=as.factor(order)), stat="boxplot",
               alpha=0.5, width=0.4, position=position_dodge(0.7)) +
  theme_bw() +
  theme(text=element_text(family="sans", face="plain", color="#000000",
                          size=15, hjust=0.5, vjust=0.5),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  guides(fill=guide_legend(title="Order")) +
  xlab("Type of Rearrangement") + ylab("Rate")
