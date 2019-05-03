source("functions.R")
library(igraph)
library(ggplot2)
library(gridExtra)
library(fdatest)
library('gprege')
library("networkD3")
library(DescTools)

#file <- "Dati/3437.edges.txt"
file <- "Dati/rip_348.edges.txt"
file<-"Dati/email-Eu-core.txt" 
#file <- "Dati/edgelistQ08.txt"
file<-"Dati/data/football.gml"

##CREATE GRAPHS
graph <- prepGraph(file,file.format="edgelist",number=TRUE) 
graph <- prepGraph(file=football,file.format="igraph") 
graph <- prepGraph(file,file.format="gml") 
graph

#Graph Random
graphRandom <- random(graph)
graphRandom


##MODULARITY
net <- cluster_fast_greedy(graph)
modularity(net)
net2 <- cluster_fast_greedy(graphRandom)
modularity(net2)

###COMMUNITY
methodCommunity (graph=graph, method="fastGreedy")#with all the method implemented in igraph
membri<-membershipCommunities(graph=graph, method="fastGreedy")
membri<-as.vector(membri)
write.csv(membri,file="Membri_email-Eu-core")

## PLOT GRAPH
plotGraph (graph)


##PLOT COMMUNITIES
plotCommu(graph,method="fastGreedy") #plot 3 D

#plot non 3 D
#fastgreedy
fg <- cluster_fast_greedy(graph)
members <- membership(fg)
V(graph)$color <- members+1
graph <- set_graph_attr(graph, "layout", layout_with_kk(graph))
plot(graph, vertex.label.dist=1.5)


##REAL RANDOM
Proc <-robinProc(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="independent")

Proc <-robinProc(graph=graph,graphRandom=graphRandom, method="edgeBetweenness",
           type="dependent")

write.csv(Proc$viMean,file="VIMeanindependent08.csv")
write.csv(Proc$ratios,file="VIratiosindependent08.csv")
write.csv(Proc$viMeanRandom,file="VIMeanRandomindependent08.csv")
write.csv(Proc$vi,file="VIindependent08.csv")
write.csv(Proc$viRandom,file="VIRandomindependent08.csv")

plotRobin(graph=graph,model=Proc$viMean,
          modelR=Proc$viMeanRandom, 
          legend=c("real data", "null model"))


##COMPARISON
Comp <- comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                method2="louvain",type="independent")

Comp <- comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                   method2="louvain",type="dependent")


plotRobinCompare(graph=graph, model1=Comp$viMean1, model2=Comp$viMean2,
                 modelR1=Comp$viMeanRandom1, modelR2=Comp$viMeanRandom2,
                 legend=c("real data", "null model"),
                 legend1vs2=c("fast greedy", "walktrap"),
                 title1="Fast Greedy",title2="walktrap",
                 title1vs2="Fast Greedy vs Walktrap")


######### TEST
#cofronto tra modello e modello nullo
robinGPTest(ratio=Proc$ratios)

robinFDATest(graph=graph, model1=Proc$viMean,model2=Proc$viMeanRandom)

robinAUCTest(graph=graph,model1=Proc$viMean,model2=Proc$viMeanRandom)


#confronto tra due metodi
robinGPTest(ratio=Comp$ratios1vs2)

robinFDATest(graph=graph, model1=Comp$viMean1, model2=Comp$viMean2, 
             legend=c("model1", "model2"))

robinAUCTest(graph=graph, model1=Comp$viMean1, model2=Comp$viMean2)


##############    PROVE   ######################
###network direzionata
###network con nomi
###network pesata
library(igraphdata)


###1
###
graph <- make_graph("Zachary")
graph <- igraph::simplify(graph) 
graphRandom <- random(graph)
List<-iter(graph=graph, graphRandom=graphRandom, method="fastGreedy",
           type="independent")
List<-iter(graph=graph, graphRandom=graphRandom, method="fastGreedy",
           type="dependent")
plotRobin(graph=graph)
Comp<-comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                 method2="walktrap",type="dependent")
Comp<-comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                 method2="walktrap",type="independent")
plotRobinCompare(graph)
karate <- make_graph("Zachary")
wc <- cluster_walktrap(karate)
members <- membership(wc)# Convert to object suitable for networkD3
karate_d3 <- igraph_to_networkD3(karate, group = members)# Create force directed network plot
forceNetwork(Links = karate_d3$links, Nodes = karate_d3$nodes,Source ='source', 
             Target ='target', NodeID ='name',Group ='group')

###2
###
data(USairports)
graph <- USairports
graph <- igraph::simplify(graph) 
graphRandom <- random(graph)
List <- iter(graph=graph,graphRandom=graphRandom, method="edgeBetweenness",
           type="independent",directed=TRUE)
#Non gira questo gira solo il comando ma con errore
# Error: Modularity is implemented for undirected graph only 
#anche se metto modularity= FALSE
graph <- prepNet(net,file.format="edgelist",method="igraph",is.directed = TRUE)
cluster_edge_betweenness(graph=graph, directed=TRUE, modularity=TRUE)


###3
###
data(Koenigsberg)
graph <- Koenigsberg
graph <- igraph::simplify(graph)
graphRandom <- random(graph)
List<-iter(graph=graph, graphRandom=graphRandom, method="fastGreedy",
           type="dependent")
#troppi pochi nodi non riesce a farlo dovremmo scriverlo 


###4
###
data(karate)
graph <- karate
graph <- igraph::simplify(graph)
graphRandom <- random(graph)
List <- iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="independent")
#non funziona 
