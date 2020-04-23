files <- list.files(path="../results/")

results.wop <- results.wp <- as.data.frame(matrix(NA, 0, 3))

colnames(results.wop) <- colnames(results.wp) <-
  c("rate", "type", "order")

for(j in 1:length(files)){
  print(files[j])
  load(file=paste("../results/",files[j],sep=""))
  if(ncol(x[[1]]) == 4){
    rates <- orders <- types <- c()
    for(i in 1:100){
      rates <- c(rates, x[[i]][26:50,2], x[[i]][26:50,3])
      types <- c(types, rep(c("asc", "desc"), each = 25))
      foo <- strsplit(files[j], split=".", fixed=T)[[1]][2]
      orders <- c(orders, rep(foo, 50))
    }
    foo <- data.frame(rates, types, orders)
    colnames(foo) <- c("rate", "type", "order")
    results.wop <- rbind(results.wop, foo)
  }
  if(ncol(x[[1]]) == 5){
    rates <- orders <- types <- c()
    for(i in 1:100){
      rates <- c(rates, x[[i]][26:50,2], x[[i]][26:50,3], x[[i]][26:50,4])
      types <- c(types, rep(c("asc", "desc", "pol"), each = 25))
      foo <- strsplit(files[j], split=".", fixed=T)[[1]][2]
      orders <- c(orders, rep(foo, 75))
    }
    foo <- data.frame(rates, types, orders)
    colnames(foo) <- c("rate", "type", "order")
    results.wp <- rbind(results.wp, foo)
  }
}
rm(list=ls()[-c(7:8)])
results.wop <- as.data.frame(results.wop)
results.wp <- as.data.frame(results.wp)

