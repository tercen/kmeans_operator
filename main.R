library(tercen)
library(dplyr)
library(reshape2)

#library(tim)
#options("tercen.workflowId" = "6015a4dd34cef273755e1a1b1500427b")
#options("tercen.stepId"     = "d31241f6-173f-473a-9307-2b4b3c5c0882")

ctx <- tercenCtx()

data <- ctx %>% 
  select(.ci, .ri, .y) %>% 
  reshape2::acast(.ci ~ .ri, value.var='.y', fill=NaN, fun.aggregate=mean) 

colnames(data) <- paste('c', colnames(data), sep='')

seed <- NULL
if(!ctx$op.value('seed') == "NULL") seed <- as.integer(ctx$op.value('seed'))


centers<-ifelse(is.null(ctx$op.value('centers')), 10, as.integer(ctx$op.value('centers')))
iter.max<-ifelse(is.null(ctx$op.value('iter.max')), 10, as.integer(ctx$op.value('iter.max')))
nstart<-ifelse(is.null(ctx$op.value('nstart')), 1, as.integer(ctx$op.value('nstart')))

#seed<-42
#centers<-2
set.seed(seed)

dataK <- kmeans(data, centers = centers , iter.max = iter.max, nstart = nstart)

df.out<-data.frame(.ci = seq(from=0,to=length(dataK$cluster)-1),
           cluster=paste0("cluster", dataK$cluster)) %>%
  ctx$addNamespace() 

df.out%>%
  ctx$save()

#tim::build_test_data(res_table = df.out, ctx = ctx, test_name = "test1")
