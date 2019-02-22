source("functions.R")
library(igraph)
library(ggplot2)
library(gridExtra)
library(fdatest)
library('gprege')
library("networkD3")
library(DescTools)

net <- "Dati/3437.edges.txt"
net <- "rip_348.edges.txt"

##CREATE GRAPHS
graph <- prepNet(net, file.format="edgelist", method="robin")
#metodo igraph un vertice in più
graphRandom <- random(graph)

##PLOT GRAPHS
plotNet(graph)

##MODULARITY
inf <- cluster_infomap(graph)
modularity(inf)
infR <- cluster_infomap(graphRandom)
modularity(infR)


fast <- cluster_fast_greedy(graph)
modularity(fast)
fastR <- cluster_fast_greedy(graphRandom)
modularity(fastR)

##PLOT COMMUNITIES
plotCommu(graph,method="infomap") #plot 3 D

#plot non 3 D
#fastgreedy
fg <- cluster_fast_greedy(graph)
members <- membership(fg)
V(graph)$color <- members+1
graph <- set_graph_attr(graph, "layout", layout_with_kk(graph))
plot(graph, vertex.label.dist=1.5)


##REAL RANDOM
List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="independent")
List<-iter(graph=graph,graphRandom=graphRandom, method="fastGreedy",
           type="dependent")

plotRobin(graph=graph)


##COMPARISON
Comp <- comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                method2="walktrap",type="independent")

Comp <- comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                method2="louvain",type="dependent")

Comp <- comparison(graph=graph,graphRandom=graphRandom,method1="fastGreedy",
                   method2="leadingEigen",type="dependent")

Comp <- comparison(graph=graph,graphRandom=graphRandom,method1="louvain",
                   method2="walktrap",type="dependent")

plotRobinCompare(graph,legend=c("real data", "null model"),
                 legend1vs2=c("fast greedy", "louvain"),
                 title1="Fast Greedy",title2="Louvain",
                 title1vs2="Fast Greedy vs Louvain")


##TEST
#cofronto tra modello e modello nullo
robinTest(graph=graph, model1=List$viMean,model2=List$viMeanRandom, 
          ratio=List$ratios, legend=c("real data", "null model"))
#confronto tra due metodi
robinTest(graph=graph, model1=Comp$viMean1,model2=Comp$viMean2, 
          ratio=Comp$ratios1vs2,legend=c("model1", "model2"))




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
