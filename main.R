source("functions.R")
library(igraph)

net <- "rip_348.edges.txt"

graph<- prepNet(net)

graphRandom <- random(graph)

iter(base="facebook_348_fastgreedy", graph, graphRandom, methodCommunity(type="fastGreedy"))

