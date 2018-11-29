source("functions.R")
library(igraph)
library(ggplot2)
net <- "rip_348.edges.txt"

graph <- prepNet(net)

graphRandom <- random(graph)

#graph<-read_graph(file="rip_348.edges.txt",format="edgelist",directed=FALSE)
#graph non viene uguale a quello di prima 

List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy", type="independent")

List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy", type="dependent")

#method<-c("walktrap", "edgeBetweenness", "fastGreedy","leadingEigen","louvain","spinglass","labelProp","infomap")

Comp<-comparison(graph=graph,method1="louvain",method2="walktrap",type="independent")



