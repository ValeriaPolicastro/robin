source("functions.R")
library(igraph)
library(ggplot2)
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
Comp<-comparison(graph=graph,method1="walktrap",method2="fastGreedy",
                 type="dependent")
Comp<-comparison(graph=graph,method1="fastGreedy",method2="walktrap",
                 type="dependent")

plotRobin(graph=graph,model1=Comp$viMean1,model2=Comp$viMean2,
          legend=c("model1", "model2"))


