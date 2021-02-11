library(tercen)
library(dplyr)
library(reshape2)

data = (ctx = tercenCtx())  %>% 
  select(.ci, .ri, .y) %>% 
  reshape2::acast(.ci ~ .ri, value.var='.y', fill=NaN, fun.aggregate=mean) 

colnames(data) = paste('c', colnames(data), sep='')


dataX = class::kmeans(data, centers=as.integer(ctx$op.value('number of clusters')))

data.frame(.ci = seq(from=0,to=length(dataX$cluster)-1),
           cluster=paste0("cluster",dataX$cluster)) %>%
  ctx$addNamespace() %>%
  ctx$save()