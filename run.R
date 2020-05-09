source("cluster_fun.R")


## read data###
data=as.matrix(read.csv("data2.csv",header = F))
pos=as.matrix(read.csv("pos.csv",header=F))
pos=as.vector(pos)

## call the cluster function
out=cluster_fun(data,pos)

clusters=out[[1]]

### each element in clusters contains row numbers of the functions in the same group

clusters



###plot the clusters
plot(out[[2]])
