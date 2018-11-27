source("functions.R")
library(igraph)
library(ggplot2)
net <- "rip_348.edges.txt"

graph <- prepNet(net)

graphRandom <- random(graph)


List<-iter(graph=graph,graphRandom=graphRandom,
           method="fastGreedy",
           type="independent",
           directed=FALSE,
           weights=NULL, 
           steps=4, 
           spins=25, 
           e.weights=NULL, 
           v.weights=NULL, 
           nb.trials=10)

List<-iter(graph=graph,graphRandom=graphRandom,method="fastGreedy", type="dependent")

#method<-c("walktrap", "edgeBetweenness", "fastGreedy","leadingEigen","louvain","spinglass","labelProp","infomap")

Comp<-comparison(graph=graph,method1="louvain",method2="walktrap",type="independent",directed=FALSE)



