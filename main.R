source("functions.R")
library(igraph)

f_net <- "rip_348.edges.txt"

g2 <- prep_net(f_net)

g2_random <- random(g2)

iter(base="facebook_348_fastgreedy",g2,g2_random,method="fastgreedy")

#modifica inutile