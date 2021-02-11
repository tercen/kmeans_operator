library(tercen)
library(dplyr)
library(reshape2)

data = (ctx = tercenCtx())  %>% 
  select(.ci, .ri, .y) %>% 
  reshape2::acast(.ci ~ .ri, value.var='.y', fill=NaN, fun.aggregate=mean) 

colnames(data) = paste('c', colnames(data), sep='')

seed <- NULL
if(!ctx$op.value('seed') == "NULL") seed <- as.integer(ctx$op.value('seed'))

dataK = kmeans(data, centers = as.integer(ctx$op.value('centers')), iter.max = as.integer(ctx$op.value('iter.max')), nstart = as.integer(ctx$op.value('nstart')))

data.frame(.ci = seq(from=0,to=length(dataK$cluster)-1),
           cluster=paste0("cluster", dataK$cluster)) %>%
  ctx$addNamespace() %>%
  ctx$save()