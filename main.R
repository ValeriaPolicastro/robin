source("functions.R")
library(igraph)

net <- "rip_348.edges.txt"

graph <- prepNet(net)

graphRandom <- random(graph)

<<<<<<< HEAD
List<-iter(base="independent",graph=graph,graphRandom=graphRandom,method="fastGreedy",type="independent")
List<-iter(base="dependent",graph=graph,graphRandom=graphRandom,method="fastGreedy",type="dependent")

plotRobin(base="dependent",graph=graph,robinList=List)

comparison(graph=graph,method1="fastGreedy",method2="louvain")
=======
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
>>>>>>> f021e7daadc6ca99c67f132cf4ac0613ed20a6d8

# writeListAsTables <- function(list, path, prefix=NULL)
# {
#     names <- names(list)
#     
#     for(name in names)
#     {
#         filaname <- file.path(paste0())
#         write.
#     }
<<<<<<< HEAD
# }
# 
#     # fileoutbats <- paste(base, "_BATS.txt", sep="")
#  filepdf <- paste(base, "_VI.pdf", sep="")
#  fileoutvirand <- paste(base, "_VI_random.txt", sep="")
#  fileoutvicase <- paste(base, "_VI_case.txt", sep="")
#  #file utilizzati da Italia per FAD
#  fileoutvirandbio <- paste(base, "_VI_random_bio.txt", sep="")
#  fileoutvicasebio <- paste(base, "_VI_case_bio.txt", sep="")


#  write.table(viRandom, fileoutvirand, sep="\t", row.names=FALSE, quote=FALSE)
#  write.table(viBhl, fileoutvicase, sep="\t", row.names=FALSE, quote=FALSE)
#  write.table(viMeanRandom, fileoutvirandbio, sep="\t", row.names=FALSE, quote=FALSE)
#  write.table(viMeanBhl, fileoutvicasebio, sep="\t", row.names=FALSE, quote=FALSE)
# #write.table(resBats, fileoutbats, sep="\t", row.names=FALSE, quote=FALSE)
# 
=======
# }
>>>>>>> f021e7daadc6ca99c67f132cf4ac0613ed20a6d8
