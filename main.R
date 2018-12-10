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
Comp<-comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                 method2="walktrap",type="dependent")

plotRobinCompare(graph)




##############    PROVE   ######################
###network direzionata
###network con nomi
###network pesata
library(igraphdata)
###1
graph<-make_graph("Zachary")
graph <- igraph::simplify(graph) 
graphRandom <- random(graph)
List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="independent")
List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="dependent")
plotRobin(graph=graph)
Comp<-comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                 method2="walktrap",type="dependent")
plotRobinCompare(graph)

###2
###
data(USairports)
graph<-USairports
graph <- igraph::simplify(graph) 
graphRandom <- random(graph)
List<-iter(graph=graph,graphRandom=graphRandom, method="edgeBetweenness",
           type="independent",directed=TRUE)
#Error: Modularity is implemented for undirected graph only

###3
###
data(Koenigsberg)
graph<-Koenigsberg
graph <- igraph::simplify(graph)
graphRandom <- random(graph)
List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="independent")
#non lo fa perchè sono nomi 



