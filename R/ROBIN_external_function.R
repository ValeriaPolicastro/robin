############### COMPARE METHOD ##########

#' robinCompare
#'
#' @description  This function compares the robustness of two community 
#' detection algorithms.
#' @param graph The output of prepGraph.
#' @param method1 The first clustering method, one of "walktrap", 
#' "edgeBetweenness", "fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","leiden","optimal","other".
#' @param args1 A \code{list} of arguments to be passed to the \code{method1} 
#' (see i.e. \link[igraph]{cluster_leiden} for a list of possible method parameters).
#' @param method2 The second custering method one of "walktrap",
#' "edgeBetweenness","fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","leiden","optimal","other".
#' @param args2 A \code{list} of arguments to be passed to the \code{method2}
#' (see i.e. \link[igraph]{cluster_leiden} for a list of possible method parameters).
#' @param FUN1 personal designed function when \code{method1} is "other". 
#' see \code{\link{methodCommunity}}.
#' @param FUN2 personal designed function when \code{method2} is "other". 
#' see \code{\link{methodCommunity}}.
#' @param measure The stability measure, one of "vi", "nmi", "split.join", 
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param ncores number of CPU cores to use.(default is 2) For a faster 
#' execution we suggest to use ncores=(parallel::detectCores(logical = FALSE)-1)
#' maximum 12 cores
#' @param type Character indicating "independent" or "dependent" for the old 
#' robin type contruction. If NULL the new faster version is computed 
#' (default NULL).
#' @param verbose flag for verbose output (default as TRUE).
#' 
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method 
#' - the matrix "Mean2" with the means of the procedure for the second method
#' 
#' @import igraph
#' @export 

#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' robinCompare(graph=graph, method1="louvain", args1 = list(resolution=0.8),
#'             method2="leiden", args2=list(objective_function ="modularity"),
#'             ncores=2)
robinCompare <-  function(graph, 
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
                          ncores=2,
                          type=NULL,
                          verbose=TRUE, distrib="NegBinom")
{
    
    methods <- c(method1, method2)
    
    # Aggiungere la versione Weigthed
    if ( is.weighted(graph) )
    {
        output <- robinCompareFastWeight(graph=graph, method1=method1, args1=args1, 
            method2=method2, args2=args2, FUN1=FUN1, FUN2=FUN2, measure=measure, 
            ncores=ncores, verbose=verbose, distrib=distrib)
    } else {
        output <- robinCompareFast(graph=graph, method1=method1, args1=args1, 
                                method2=method2, args2=args2, 
                                FUN1=FUN1, FUN2=FUN2, measure=measure, 
                                ncores=ncores, verbose=verbose)
    }
        
     if(any(type %in% c("independent", "dependent")))
     {
    
          output <- robinCompareNoParallel(graph=graph, method1=method1, args1=args1,
                            method2=method2, args2=args2, measure=measure, 
                           type=type) 
     }
   
   
    outputRobin <- c(output, model=methods, list(graph=graph))

   class(outputRobin) <- "robin"
    return(outputRobin)
    
}


############### ROBUST METHOD ##########
#' robinRobust
#' @description This functions implements a procedure to examine the stability 
#' of the partition recovered by some algorithm against random perturbations 
#' of the original graph structure.
#' @param graph The output of prepGraph.
#' @param graphRandom The output of random function.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "leiden","optimal".
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param measure The stability measure, one of "vi", "nmi", "split.join", 
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param type The type of robin construction, dependent or independent 
#' procedure.
#' @param ... other parameter.
#' @param verbose flag for verbose output (default as TRUE).
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean" with the means of the procedure for the graph
#' - the matrix "MeanRandom" with the means of the procedure for the random graph. 
#' 
#' @import igraph
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' robinRobust(graph=graph, graphRandom=graphRandom, method="louvain",
#'     resolution=0.8, measure="vi", type="independent")


robinRobust <-  function(graph, graphRandom, 
                          method=c("walktrap", "edgeBetweenness", 
                                   "fastGreedy", "louvain", "spinglass", 
                                   "leadingEigen", "labelProp", "infomap",
                                   "optimal", "leiden", "other"),
                          ...,
                          FUN=NULL, measure= c("vi", "nmi","split.join", "adjusted.rand"),
                          type=c("independent","dependent"), nrep=5,verbose=TRUE )
{
    type <- match.arg(type)
    methods <- c("real data", "null model")
     # No Parallel              
    output <- robinRobustNoParallel(graph=graph, graphRandom= graphRandom, 
                                    method=method,
                                    ...,
                                    FUN=NULL, measure=measure,
                                    type=type, nrep=nrep, verbose=verbose)
     # Parallel version: 
    #DA FARE

outputRobin <- c(output, model=methods, list(graph=graph))

class(outputRobin) <- "robin"
return(outputRobin)
}
                          

                          
            

