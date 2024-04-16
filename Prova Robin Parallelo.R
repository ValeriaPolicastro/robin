proc1F <- robinRobustFast(graph=graph, graphRandom=graphRandom, method="louvain",
                          resolution = 0.8, measure="vi")
proc2F <- robinRobustFast(graph=graph, graphRandom=graphRandom, method="leiden",
                          objective_function = "modularity", measure="vi")


proc3F <- robinRobustFast(graph=graph, graphRandom=graphRandom, method="infomap",
                          modularity = FALSE, measure="vi")



# Infomap e Leiden fast vanno in errore

proc1 <- robinRobustNoParallel(graph=graph, graphRandom=graphRandom, method="louvain",
                               resolution = 0.8, measure="vi", type="independent")

proc2 <- robinRobustNoParallel(graph=graph, graphRandom=graphRandom, method="leiden",
                               objective_function = "modularity", measure="vi", type="independent")


proc3 <- robinRobustNoParallel(graph=graph, graphRandom=graphRandom, method="infomap",
                               modularity = FALSE, measure="vi", type="independent")
