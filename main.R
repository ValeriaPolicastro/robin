source("functions.R")
library(igraph)
library(ggplot2)
library(gridExtra)
net <- "rip_348.edges.txt"

##CREATE GRAPHS
graph <- prepNet(net,file.format="edgelist",method="robin")
#metodo igraph un vertice in più
graphRandom <- random(graph)


##REAL RANDOM
List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="independent")
List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="dependent")

plotRobin(graph=graph)


##COMPARISON
Comp<-comparison(graph=graph,graphRandom=graphRandom,method1="walktrap",
                 method2="fastGreedy",type="independent")
Comp<-comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",method2="walktrap",
                 type="dependent")

plotRobinCompare(graph,legend1vs2=c("method1", "method2"))

