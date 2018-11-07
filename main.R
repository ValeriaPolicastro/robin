source("functions.R")
library(igraph)

net <- "rip_348.edges.txt"

graph<- prepNet(net)

graphRandom <- random(graph)

iter(base="output",graph=graph,graphRandom=graphRandom,method="fastGreedy",type="independent")

comparison(method1,method2)

