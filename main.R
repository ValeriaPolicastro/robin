source("functions.R")
library(igraph)

net <- "rip_348.edges.txt"

graph <- prepNet(net)

graphRandom <- random(graph)

iterList <- iter(base="output",
    graph=graph,
    graphRandom=graphRandom,
    method="fastGreedy", 
    type="independent")



iterList <- iter(base="output",
                 graph=graph,
                 graphRandom=graphRandom,
                 method="fastGreedy", 
                 type="dependent")



comparison(method1, method2)

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