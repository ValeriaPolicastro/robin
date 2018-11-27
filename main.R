source("functions.R")
library(igraph)
library(ggplot2)
net <- "rip_348.edges.txt"

graph <- prepNet(net)

graphRandom <- random(graph)


List<-iter(graph=graph,graphRandom=graphRandom,method="edgeBetweenness",type="independent")

List<-iter(graph=graph,graphRandom=graphRandom,method="edgeBetweenness",type="dependent")

#method<-c("walktrap", "edgeBetweenness", "fastGreedy","leadingEigen","louvain","spinglass","labelProp","infomap")

Comp<-comparison(graph=graph,method1="louvain",method2="walktrap",type="independent")




# writeListAsTables <- function(list, path, prefix=NULL)
# {
#     names <- names(list)
#     
#     for(name in names)
#     {
#         filaname <- file.path(paste0())
#         write.
#     }

# }
# 
