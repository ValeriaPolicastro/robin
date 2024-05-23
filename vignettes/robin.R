## -----------------------------------------------------------------------------
#install.packages("robin")

## -----------------------------------------------------------------------------
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("gprege")
# 
# install.packages("robin")

## ----message=FALSE, warning=FALSE, paged.print=TRUE---------------------------
library("robin")

## -----------------------------------------------------------------------------
my_network <- system.file("example/football.gml", package="robin")
# downloaded from: http://www-personal.umich.edu/~mejn/netdata/
graph <- prepGraph(file=my_network, file.format="gml")
graph

## -----------------------------------------------------------------------------
plotGraph(graph)

## -----------------------------------------------------------------------------
methodCommunity(graph=graph, method="fastGreedy") 

## -----------------------------------------------------------------------------
membershipCommunities(graph=graph, method="fastGreedy") 

## -----------------------------------------------------------------------------
members <- membershipCommunities(graph=graph, method="fastGreedy")
plotComm(graph=graph, members=members)

## -----------------------------------------------------------------------------
graphRandom <- random(graph=graph)
graphRandom

