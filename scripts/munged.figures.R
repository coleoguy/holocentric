library(ggraptR)
library(ggplot2)
#ggraptR(results.wop)

#Plot of results without polyploidy
ggplot(results.wop, aes(y=rate, x=as.factor(type))) + 
  geom_boxplot(aes(fill=as.factor(order)), stat="boxplot", 
               alpha=0.5, width=0.2, position=position_dodge(0.5)) + 
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", 
                          size=15, hjust=0.5, vjust=0.5), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  guides(fill=guide_legend(title="Order")) + 
   xlab("Type of Rearrangement") + ylab("Rate")

#Plot of results with polyploidy
ggplot(results.wp, aes(y=rate, x=as.factor(type))) + 
  geom_boxplot(aes(fill=as.factor(order)), stat="boxplot", 
               alpha=0.5, width=0.2, position=position_dodge(0.5)) + 
  theme_bw() + 
  theme(text=element_text(family="sans", face="plain", color="#000000", 
                          size=15, hjust=0.5, vjust=0.5), 
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  guides(fill=guide_legend(title="Order")) + 
  xlab("Type of Rearrangement") + ylab("Rate")
