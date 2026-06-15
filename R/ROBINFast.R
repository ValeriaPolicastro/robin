#' robinCompareFast
#'
#' @description This function compares two community detection algorithms.
#' Is the parallelized and faster version of \code{\link{robinCompare}}
#'
#' @param graph The output of prepGraph.
#' @param method1 The first clustering method, one of "walktrap", 
#' "edgeBetweenness", "fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal".
#' @param args1 A \code{list} of arguments to be passed to the \code{method1} 
#' (see i.e. \link[igraph]{cluster_leiden} for a list of possible method parameters).
#' @param method2 The second custering method one of "walktrap",
#' "edgeBetweenness","fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal".
#' @param args2 A \code{list} of arguments to be passed to the \code{method2}
#' (see i.e. \link[igraph]{cluster_leiden} for a list of possible method parameters).
#' @param FUN1 personal designed function when \code{method1} is "others". 
#' see \code{\link{methodCommunity}}.
#' @param FUN2 personal designed function when \code{method2} is "others". 
#' see \code{\link{methodCommunity}}.
#' @param measure The stability measure, one of "vi", "nmi", "split.join", 
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param verbose flag for verbose output (default as TRUE).
#' @param seed set seed (default seed=NULL)
#' @param BPPARAM optional BiocParallel parameter object. If \code{NULL}
#' (default), uses \code{parallel::mclapply} on Unix/macOS or sequential
#' \code{lapply} on Windows. If a BiocParallel object is supplied and the
#' \pkg{BiocParallel} package is installed, it will be used instead.
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method 
#' - the matrix "Mean2" with the means of the procedure for the second method
#' 
#' @import igraph
#' @keywords internal
# @examples my_file <- system.file("example/football.gml", package="robin")
# graph <- prepGraph(file=my_file, file.format="gml")
# robinCompareFast(graph=graph, method1="louvain", args1 = list(resolution=0.8),
# method2="leiden", args2=list(leiden_objective_function ="modularity"))

robinCompareFast <- function(graph, 
                         method1=c("walktrap", "edgeBetweenness", "fastGreedy",
                                   "leadingEigen","louvain","spinglass",
                                   "labelProp","infomap","optimal","leiden",
                                   "other"),
                         args1=list(),
                         method2=c("walktrap", "edgeBetweenness", "fastGreedy",
                                   "leadingEigen","louvain","spinglass",
                                   "labelProp","infomap","optimal","leiden",
                                   "other"),
                         args2=list(),
                         measure= c("vi", "nmi", "split.join", "adjusted.rand"),
                         FUN1=NULL, FUN2=NULL, seed=NULL,
                         verbose=TRUE, BPPARAM=NULL)
{   
    method1 <- match.arg(method1)
    method2 <- match.arg(method2)
    measure <- match.arg(measure)
    args11 <- c(list(graph=graph), method=method1, FUN=FUN1, args1)
    args21 <- c(list(graph=graph), method=method2, FUN=FUN2, args2)
    if(!is.null(seed)){set.seed(seed)}
    comReal1 <- do.call(robin::membershipCommunities, args11)
    N <- igraph::vcount(graph)
    
    n_communities1 <- length(table(comReal1))
    
    if (n_communities1 == N) {
        warning("Each community has only one node with method1")
    }
    
    if (n_communities1 == 1) {
        warning("Only one community with method1")
    }
    
    comReal2 <- do.call(robin::membershipCommunities, args21)
   
    n_communities2 <- length(table(comReal2))
    
    if (n_communities2 == N) {
        warning("Each community has only one node with method2")
    }
    
    if (n_communities2 == 1) {
        warning("Only one community with method2")
    }
    
   
    de <- igraph::gsize(graph)
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detected robin method type independent\nIt can take time ... It depends on the size of the network.\n")
    vet1 <- seq(5, 60, 5) 
    vet <- round(vet1*de/100, 0)
   
     parfunct <- function(z, graph, method1, method2, comReal1, comReal2, N, 
                          measure, args1, args2, FUN1, FUN2)
     {
        # print(list(args1,args2))
        
        
        
        MeansList <- lapply(1:10, function(s)
        {
            
            graphRList <- igraph::rewire(graph, 
                                with=igraph::keeping_degseq(loops=FALSE,
                                                             niter=z))
          
            argsP <- c(list(graph=graphRList), method=method1, FUN=FUN1, args1)
            comr1 <- do.call(robin::membershipCommunities, argsP)
          
            
            
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
                                               method=measure)))/2
                
            }else{
                measure1 <- 1-(igraph::compare(comr1, comReal1, 
                                               method=measure))
                
            }
            
            
            argsP <- c(list(graph=graphRList), method=method2, FUN=FUN2, args2)
            comr2 <- do.call(robin::membershipCommunities, argsP)
          
            if(measure=="vi")
            {
                
                measure2 <-  igraph::compare(comr2, comReal2, 
                                             method=measure)/log2(N)
            } else if (measure=="split.join"){
                measure2 <- igraph::compare(comr2, comReal2, 
                                            method=measure)/(2*N)
            } else if (measure=="adjusted.rand"){
                
                measure2 <- (1-(igraph::compare(comr2, comReal2, 
                                               method=measure)))/2
            }else{
                
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
    
    
    zlist <- robin_lapply(vet, parfunct, graph=graph, measure=measure,
                          method1=method1, method2=method2, args1=args1, args2=args2,
                          comReal1=comReal1, N=N, comReal2=comReal2, FUN1=FUN1,
                          FUN2=FUN2, BPPARAM=BPPARAM)
    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))
    
    Measure1 <- cbind(rep(0, 10), Measure1)
    Measure2 <- cbind(rep(0, 10), Measure2)
    
    colnames(Measure1) <- nRewire 
    colnames(Measure2) <- nRewire 
    return(list(Mean1=Measure1,
                Mean2=Measure2,
                Communities1=comReal1,
                Communities2=comReal2))
}




#' robinRobustFast
#'
#' @description This functions implements a procedure to examine the stability 
#' of the partition recovered by some algorithm against random perturbations 
#' of the original graph structure.
#'
#' @param graph The output of prepGraph.
#' @param graphRandom The output of random function.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "leiden","optimal".
#' @param measure The stability measure, one of "vi", "nmi", "split.join", 
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param FUN1 in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param ... other parameter
#' @param verbose flag for verbose output (default as TRUE)
#' @param seed set seed (default seed=NULL)
#' @param BPPARAM optional BiocParallel parameter object. If \code{NULL}
#' (default), uses \code{parallel::mclapply} on Unix/macOS or sequential
#' \code{lapply} on Windows. If a BiocParallel object is supplied and the
#' \pkg{BiocParallel} package is installed, it will be used instead.
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean" with the means of the procedure for the graph
#' - the matrix "MeanRandom" with the means of the procedure for the random graph. 
#' 
#' @import igraph
#' @keywords internal
# @example my_file <- system.file("example/football.gml", package="robin")
# graph <- prepGraph(file=my_file, file.format="gml")
# graphRandom <- random(graph=graph)
# robinRobustFast(graph=graph, graphRandom=graphRandom, method="leiden",
# resolution_parameter = 1, measure="vi")
robinRobustFast <- function(graph, graphRandom, 
                            method=c("walktrap", "edgeBetweenness", 
                                     "fastGreedy", "louvain", "spinglass", 
                                     "leadingEigen", "labelProp", "infomap",
                                     "optimal", "leiden", "other"), ..., FUN1=NULL, 
                            measure=c("vi", "nmi", "split.join","adjusted.rand"),
                            seed=NULL,
                            verbose=TRUE, BPPARAM=NULL)
{
    method <- match.arg(method)
    measure <- match.arg(measure)
    N <- igraph::vcount(graph)
    if(!is.null(seed)){set.seed(seed)}
    comReal1 <- membershipCommunities(graph=graph, method=method,
                                      FUN=FUN1, ...=...)
    
     n_communitiesReal <- length(table(comReal1))
    
    if (n_communitiesReal == N) {
        warning("Each community has only one node in the graph")
    }
    
    if (n_communitiesReal == 1) {
        warning("Only one community in the graph")
    }
    
    comReal2 <- membershipCommunities(graph=graphRandom, method=method,
                                       FUN=FUN1, ...=...)
    # random network
    n_communitiesRandom <- length(table(comReal2))
    
    if (n_communitiesRandom == N) {
        warning("Each community has only one node in the graphRandom")
    }
    
    if (n_communitiesRandom == 1) {
        warning("Only one community in the graphRandom")
    }
    
   
    de <- igraph::gsize(graph)
    nRewire <- seq(0, 60, 5)
    if(verbose) cat("Detected robin method type independent\nIt can take time ... It depends on the size of the network.\n")
    vet1 <- seq(5, 60, 5) 
    vet <- round(vet1*de/100, 0)
    
    parfunct <- function(z, graph, method, comReal1, comReal2, N, 
                         measure, FUN1, ...)
    {
        # print(list(...))
        MeansList <- lapply(1:10, function(s)
        {
            
            graphRList <- igraph::rewire(graph, 
                                         with=igraph::keeping_degseq(loops=FALSE,
                                                                     niter=z))
            
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
            }else if (measure=="adjusted.rand") {
                measure1 <- (1-(igraph::compare(comr1, comReal1, 
                                               method=measure)))/2
                
            }else{
                measure1 <- 1-(igraph::compare(comr1, comReal1, 
                                               method=measure))
            }
            
            graphRandomList <- igraph::rewire(graphRandom, 
                                         with=igraph::keeping_degseq(loops=FALSE,
                                                                     niter=z))
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
            } else if (measure=="adjusted.rand"){
                
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
    
    zlist <- robin_lapply(vet, FUN=parfunct, graph=graph,
                          method=method, comReal1=comReal1, comReal2=comReal2, N=N,
                          measure=measure, FUN1=FUN1, ...=..., BPPARAM=BPPARAM)

    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))
    
    Measure1 <- cbind(rep(0, 10), Measure1)
    Measure2 <- cbind(rep(0, 10), Measure2)
    
    colnames(Measure1) <- nRewire 
    colnames(Measure2) <- nRewire 
    return(list(Mean=Measure1,
                MeanRandom=Measure2,
                CommunitiesReal=comReal1,
                CommunitiesRandom=comReal2))
}


