####### GRAPH RANDOM WEIGHTED#########
#' randomWeight
#'
#' @description This function randomly rewires the edges while preserving the original graph's
#' degree distribution.
#' @param graph The output of prepGraph.
#' @param dist Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "Other". See
#'   \code{\link[perturbR]{rewireR}} for details.
#' @param verbose flag for verbose output (default as FALSE)
#'
#' @return An igraph object, a randomly rewired graph.
#' @import igraph
#' @keywords internal
#'

randomWeight <- function(graph, dist="Other", verbose=FALSE)
{
    if(verbose) cat("Randomizing the graph edges.\n")
    v <- igraph::vcount(graph) ## number of vertex
    numberPerturbAll <- round((v*(v-1))/2, 0)
    adj <- as_adjacency_matrix(graph, attr="weight", sparse=FALSE)
    graphRandom <- as.matrix(perturbR::rewireR(sym.matrix=adj, 
                                               nperturb=numberPerturbAll, 
                                               dist=dist))
    #rewiring for z all the edges
    graphRandom <- graph_from_adjacency_matrix(graphRandom, weighted=TRUE, 
                                               mode="undirected")
    return(graphRandom)
}


 ####### REWIRE WEIGHTED Internal #########
# 
# #' rewireWeight
# #' @description makes the rewire for weighted networks
# #' @param data The output of prepGraph
# #' @param number Number of rewiring trials to perform.
# #' @param type method to rewire weighted graphs
# #' @keywords internal
# rewireWeight <- function(data, number, type=NULL)
# {
#     if(type=="Shuffle")
#     {
#         print("Shuffle method")
#            <- sample(E(graph)$weight)
#            
#            
#          
#     }else if(type=="Garlaschelli method"){
#         print("Garlaschelli Method")
#         
#      
#              }else if(type=="Sum method"){
#        
#         print("Keep Sum and distribution weight Method") 
# 
# 
# n_numbers <- E(graph)$weight
# 
# # Subset di p numeri (p <= length(n_numbers))
# p <- number
# p_subset <- sample(n_numbers, p)
# 
# # Funzione per ottenere un nuovo subset di p numeri mantenendo la somma e la distribuzione
# generate_subset <- function(n_numbers, p_subset) {
#     target_sum <- sum(p_subset) # La somma target da mantenere
#     n_dist <- n_numbers / sum(n_numbers) # Distribuzione percentuale del set originale
# 
#     # Campiona p numeri dalla distribuzione di n_numbers
#     sampled_indices <- sample(1:length(n_numbers), p, replace = TRUE)
#     sampled_dist <- n_dist[sampled_indices]
# 
#     # Calcola i nuovi valori proporzionali alla distribuzione e alla somma target
#     new_subset <- target_sum * (sampled_dist / sum(sampled_dist))
# 
#     # Aggiusta leggermente i valori per assicurarsi che la somma sia esattamente quella target
#     diff <- target_sum - sum(new_subset)
#     adjustment <- diff / p
#     new_subset <- new_subset + adjustment
# 
#     return(new_subset)
# }
# 
# # Applica la funzione
# new_p_subset <- generate_subset(n_numbers, p_subset)
# 
# # Stampa il risultato
# print(new_p_subset)
# print(sum(new_p_subset)) # Dovrebbe essere uguale a sum(p_subset)
# 
# 
#     }else {
#     print("Rewire weight method")
#     graphRewire <- igraph::rewire(data, with=keeping_degseq(loops=FALSE,
#                                                             niter=number))
#      NotChaged <- igraph::intersection(graph, graphRewire)
#      newWeight <- sample(E(difference(graph,graphRewire))$weight)
#     EdgeAggiunti <- E(difference(graphRewire,graph))
#      gg <- difference(graphRewire,graph)
#      E(gg)$weight <- newWeight
#     U <- union(gg,NotChaged)
#      E(U)$weight_1[which(is.na(E(U)$weight_1), arr.ind = TRUE)] <- E(U)$weight[which(!is.na(E(U)$weight),
#                                                                                      arr.ind = TRUE)]
#      E(U)$weight <- E(U)$weight_1
#      U <- delete_edge_attr(U, "weight_1")
#      U <- delete_edge_attr(U, "weight_2")
#      }
#     return(U)
# }


########## ROBIN COMPARE WEIGHTED#############

#' robinCompareFastWeight
#'
#' @description This function compares two community detection algorithms, from
#' weighted networks. Is the parallelized and faster version.
#' @param graph The output of prepGraph.
#' @param method1 The first clustering method, one of "walktrap",
#' "edgeBetweenness", "fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal","leiden".
#' @param args1 A \code{list} of arguments to be passed to the \code{method1} 
#' (see i.e. \link[igraph]{cluster_leiden} for a list of possible method parameters).
#' @param method2 The second custering method one of "walktrap",
#' "edgeBetweenness","fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal","leiden".
#' @param args2 A \code{list} of arguments to be passed to the \code{method2}
#' (see i.e. \link[igraph]{cluster_leiden} for a list of possible method parameters).
#' @param measure The stability measure, one of "vi", "nmi", "split.join",
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param FUN1 personal designed function when method1 is "others".
#' see \code{\link{methodCommunity}}.
#' @param FUN2 personal designed function when method2 is "others".
#' see \code{\link{methodCommunity}}.
#' @param verbose flag for verbose output (default as TRUE).
#' @param dist Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "Other". See \link[perturbR]{rewireR}
#'  for details.
#' @param BPPARAM the BiocParallel object of class \code{bpparamClass} that 
#' specifies the back-end to be used for computations. See
#' \link[BiocParallel]{bpparam} for details.
#'
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method
#' - the matrix "Mean2" with the means of the procedure for the second method
#'
#' @import igraph parallel perturbR
#' @importFrom BiocParallel bplapply bpparam
#' @keywords internal

robinCompareFastWeight <- function(graph,
                                   method1=c("walktrap", "edgeBetweenness", "fastGreedy",
                                             "leadingEigen", "louvain", "spinglass",
                                             "labelProp", "infomap", "optimal", "leiden", 
                                             "other"),
                                   args1=list(),
                                   method2=c("walktrap", "edgeBetweenness", "fastGreedy",
                                             "leadingEigen", "louvain", "spinglass",
                                             "labelProp", "infomap", "optimal", "leiden", 
                                             "other"),
                                   args2=list(),
                                   FUN1=NULL, FUN2=NULL,
                                   measure=c("vi", "nmi","split.join", "adjusted.rand"),
                                   #ncores=2,
                                   verbose=TRUE, dist="Other", BPPARAM=BiocParallel::bpparam())
{
    method1 <- match.arg(method1)
    method2 <- match.arg(method2)
    measure <- match.arg(measure)
    args11 <- c(list(graph=graph), method=method1, FUN=FUN1, args1)
    args21 <- c(list(graph=graph), method=method2, FUN=FUN2, args2)
    comReal1 <- do.call(robin::membershipCommunities, args11)
    comReal2 <- do.call(robin::membershipCommunities, args21)
    N <- igraph::vcount(graph)
    de <- round((N*(N-1))/2, 0)
    Measure <- NULL
    vector1 <- NULL
    vector2 <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detected robin method type independent\nIt can take time ... It depends on the size of the network.\n")
    vet1 <- seq(5, 60, 5)
    vet <- round(vet1*de/100, 0)
     #print(vet)
    parfunct <- function(z, graph, method1, method2, comReal1, comReal2, N, 
                         measure, args1, args2, FUN1, FUN2,dist)
    {
        #print(list(args1,args2))
        
        #print(z)
        
        MeansList <- lapply(1:10, function(s)
        {
    
           #print(s)
            #print(z)
             adj <- igraph::as_adjacency_matrix(graph, attr="weight", sparse = FALSE)
            gR <- as.matrix(perturbR::rewireR(adj, z,dist = dist))
            graphRList <- igraph::graph_from_adjacency_matrix(gR,weighted = TRUE,
                                                              mode="undirected")
            
            argsP <- c(list(graph=graphRList), method=method1, FUN=FUN1, args1)
            comr1 <- do.call(membershipCommunities, argsP)
         
            
            if(measure=="vi")
            {
                measure1 <- igraph::compare(comr1, comReal1,
                                            method=measure)/log2(N)
            } else if(measure=="split.join")
            {
                measure1 <- igraph::compare(comr1, comReal1,
                                            method=measure)/(2*N)
            }else if(measure=="adjusted.rand"){
                measure1 <- (1-(igraph::compare(comr1, comReal1,
                                               method=measure)))/2}
            else{
                measure1 <- 1-(igraph::compare(comr1, comReal1,
                                               method=measure))
                
            }
            
            
            argsP <- c(list(graph=graphRList), method=method2, FUN=FUN2, args2)
            comr2 <- do.call(membershipCommunities, argsP)
           
            if(measure=="vi")
            {
                
                measure2 <-  igraph::compare(comr2, comReal2,
                                             method=measure)/log2(N)
            } else if (measure=="split.join"){
                measure2 <- igraph::compare(comr2, comReal2,
                                            method=measure)/(2*N)
            } else if(measure=="adjusted.rand"){
                
                measure2 <- (1-(igraph::compare(comr2, comReal2,
                                               method=measure)))/2
            } else{
                
                measure2 <- 1-(igraph::compare(comr2, comReal2,
                                               method=measure))
            }
            return(list("Measure1"=measure1, "Measure2"=measure2))
        })

        
        m1 <- unlist(lapply(MeansList, function(mm)
        {
            mm$Measure1
        }))
        
        m2 <- unlist(lapply(MeansList, function(mm)
        {
            mm$Measure2
        }))
        
    
        return(list("Measure1"=m1,"Measure2"=m2))
    }
    
    zlist <- BiocParallel::bplapply(vet, parfunct, graph=graph, measure=measure,
                                    method1=method1, method2=method2, args1=args1, args2=args2,
                                    comReal1=comReal1, N=N, comReal2=comReal2, FUN1=FUN1,
                                    FUN2=FUN2, dist=dist, BPPARAM=BPPARAM)
    
    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))

    Measure1 <- cbind(rep(0, 10), Measure1)
    Measure2 <- cbind(rep(0, 10), Measure2)
    
    colnames(Measure1) <- nRewire
    colnames(Measure2) <- nRewire
    return(list(Mean1=Measure1,
                Mean2=Measure2))
}




########### ROBIN ROBUST WEIGHTED #########

#' robinRobustFastWeighted
#'
#' @description This functions implements a procedure to examine the stability
#' of the partition recovered by some algorithm against random perturbations
#' of the original graph structure for weighted network.
#' @param graph The output of prepGraph.
#' @param graphRandom The output of random function.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness",
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "leiden","optimal".
#' @param FUN1 in case the @method parameter is "other" there is the possibility
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights
#' (that can be NULL), and has to return a community object.
#' @param measure The stability measure, one of "vi", "nmi", "split.join",
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param ... other parameter
#' @param dist Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "Other". See \link[perturbR]{rewireR}
#' for details.
#' @param verbose flag for verbose output (default as TRUE).
#' @param BPPARAM the BiocParallel object of class bpparamClass that 
#' specifies the back-end to be used for computations. See 
#' \link[BiocParallel]{bpparam} for details.
#'
#' @return A list object with two matrices:
#' - the matrix "Mean" with the means of the procedure for the graph
#' - the matrix "MeanRandom" with the means of the procedure for the random graph.
#' @keywords internal
#' @import igraph parallel perturbR
#' @importFrom BiocParallel bplapply bpparam

robinRobustFastWeighted <- function(graph, graphRandom, 
                                    method=c("walktrap", "edgeBetweenness",
                                             "fastGreedy", "louvain", "spinglass",
                                             "leadingEigen", "labelProp", "infomap",
                                             "optimal", "leiden", "other"),
                                    ..., FUN1=NULL,
                                    measure= c("vi", "nmi", "split.join", "adjusted.rand"),
                                    verbose=TRUE, dist="Other",
                                    BPPARAM=BiocParallel::bpparam())
{   
    method <- match.arg(method)
    measure <- match.arg(measure)
    comReal1 <- membershipCommunities(graph=graph, method=method,
                                      FUN=FUN1, ...=...) 
    comReal2 <- membershipCommunities(graph=graphRandom, method=method,
                                      FUN=FUN1, ...=...)
    de <- igraph::gsize(graph)
    N <- igraph::vcount(graph)
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detected robin method type independent\nIt can take time ... It depends on the size of the network.\n")
    vet1 <- seq(5, 60, 5) 
    vet <- round(vet1*de/100, 0)
    
    parfunct <- function(z, graph, method, comReal1, comReal2, N, 
                         measure, dist, FUN1, ...)
    {
        #print(list(...))
        MeansList <- lapply(1:10, function(s)
        {
            
            adj <- igraph::as_adjacency_matrix(graph, attr="weight", sparse = FALSE)
            gR <- as.matrix(perturbR::rewireR(adj, z,dist = dist))
            graphRList <- igraph::graph_from_adjacency_matrix(gR,weighted = TRUE,
                                                              mode="undirected")
            
            
            comr1 <- robin::membershipCommunities(graph=graphRList,
                                                  method=method,
                                                  FUN=FUN1,
                                                  ...=...)
            
            if(measure=="vi")
            {
                measure1 <- igraph::compare(comr1, comReal1, 
                                            method=measure)/log2(N)
            } else if(measure=="split.join")
            {
                measure1 <- igraph::compare(comr1, comReal1, 
                                            method=measure)/(2*N)
            }else if (measure=="adjusted.rand"){
                measure1 <- (1-(igraph::compare(comr1, comReal1, 
                                               method=measure)))/2}
            else{
                measure1 <- 1-(igraph::compare(comr1, comReal1, 
                                               method=measure))
                
            }
            
            adj <- igraph::as_adjacency_matrix(graphRandom, attr="weight", sparse = FALSE)
            gR <- as.matrix(perturbR::rewireR(adj, z,dist = dist))
            graphRandomList <- igraph::graph_from_adjacency_matrix(gR,weighted = TRUE,
                                                              mode="undirected")
            
          
            comr2 <- robin::membershipCommunities(graph=graphRandomList, 
                                                  FUN=FUN1,
                                                  method=method,
                                                  ...=...)
            if(measure=="vi")
            {
                
                measure2 <-  igraph::compare(comr2, comReal2, 
                                             method=measure)/log2(N)
            } else if (measure=="split.join"){
                measure2 <- igraph::compare(comr2, comReal2, 
                                            method=measure)/(2*N)
            }else if(measure=="adjusted.rand"){
                
                measure2 <- (1-(igraph::compare(comr2, comReal2, 
                                               method=measure)))/2
            } else{
                
                measure2 <- 1-(igraph::compare(comr2, comReal2, 
                                               method=measure))
            }
            return(list("Measure1"=measure1, "Measure2"=measure2))
        })
        
        
        m1 <- unlist(lapply(MeansList, function(mm)
        {
            mm$Measure1
        }))
        
        m2 <- unlist(lapply(MeansList, function(mm)
        {
            mm$Measure2
        }))
        
        
        
        
        
        return(list("Measure1"=m1,"Measure2"=m2))
    }
    zlist <- BiocParallel::bplapply(vet, parfunct, graph=graph, measure=measure,
                                    method=method, comReal1=comReal1, N=N,
                                    FUN1=FUN1, comReal2=comReal2, dist=dist,
                                    ...=...,BPPARAM=BPPARAM)
    
    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))
    
    Measure1 <- cbind(rep(0, 10), Measure1)
    Measure2 <- cbind(rep(0, 10), Measure2)
    
    colnames(Measure1) <- nRewire 
    colnames(Measure2) <- nRewire 
    return(list(Mean=Measure1,
                MeanRandom=Measure2))
    }


