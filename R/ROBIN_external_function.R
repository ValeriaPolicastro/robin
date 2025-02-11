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
#' @param type The type of robin construction, dependent or independent.
#' @param dist Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "Other". See
#'   \code{\link[perturbR]{rewireR}} for details.
#' @param BPPARAM the BiocParallel object of class \code{bpparamClass} that 
#' specifies the back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
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
#'             method2="leiden")
## Weighted example:
# E(graph)$weight <- round(runif(ecount(graph),min=1,max=10))
# robinCompare(graph=graph, method1="louvain", args1 = list(resolution=0.8), 
# method2="leiden")

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
                          type="independent",
                          verbose=TRUE,rewire.w.type=c("Rewire","Shuffle","Garlaschelli","Sum"),
                          dist="Other",BPPARAM=BiocParallel::bpparam())
{
    
    methods <- c(method1, method2)
    
    # Weigthed version
    if ( is_weighted(graph) )
    {
       print("Weighted Network Parallel Function")
         output <- robinCompareFastWeight(graph=graph, method1=method1, args1=args1, 
            method2=method2, args2=args2, FUN1=FUN1, FUN2=FUN2, measure=measure, 
            verbose=verbose, dist=dist,rewire.w.type=rewire.w.type, BPPARAM=BPPARAM)
    } else {
        if(type=="dependent")
        {
            
           print("Unweighted Network No Parallel Function")
             output <- robinCompareNoParallel(graph=graph, method1=method1, args1=args1,
                                             method2=method2, args2=args2, measure=measure, 
                                             type=type) 
        }else{
            print("Unweighted Network Parallel Function")
            output <- robinCompareFast(graph=graph, method1=method1, args1=args1, 
                                       method2=method2, args2=args2, 
                                       FUN1=FUN1, FUN2=FUN2, measure=measure, 
                                       verbose=verbose,BPPARAM=BPPARAM)
            
        }
        
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
#' @param ... other parameters for the community detection methods.
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param measure The stability measure, one of "vi", "nmi", "split.join", 
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param type The type of robin construction, dependent or independent.
#' @param rewire.w.type for weighted graph. Option to rewire one of "Rewire",
#' "Shuffle","Garlaschelli","Sum" 
#' @param dist for weighted graph with "Garlaschelli" @rewire.w.type method. 
#' Option to rewire in a manner that retains overall graph weight regardless of 
#' distribution of edge weights. This option is invoked by putting any text into
#'  this field. Defaults to "Other". See \code{\link[perturbR]{rewireR}} for 
#'  details.
#' @param BPPARAM the BiocParallel object of class \code{bpparamClass} that 
#' specifies the back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
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
#' robinRobust(graph=graph, graphRandom=graphRandom, method="leiden")
##    Weighted Example:
# E(graph)$weight <- round(runif(ecount(graph),min=1,max=10))
# graphRandom <- random(graph=graph)
# robinRobust(graph=graph, graphRandom=graphRandom, method="leiden")

robinRobust <-  function(graph, graphRandom, 
                          method=c("walktrap", "edgeBetweenness", 
                                   "fastGreedy", "louvain", "spinglass", 
                                   "leadingEigen", "labelProp", "infomap",
                                   "optimal", "leiden", "other"),
                          ...,
                          FUN=NULL, measure= c("vi", "nmi","split.join", "adjusted.rand"),
                         type="independent",verbose=TRUE,
                         rewire.w.type=c("Rewire","Shuffle","Garlaschelli","Sum"),
                         dist="Other",
                             BPPARAM=BiocParallel::bpparam())
{

    methods <- c("real data", "null model")
    
    # Weigthed version
    if ( is_weighted(graph) )
    {
        print("Weighted Network Parallel Function")
        output <- robinRobustFastWeighted(graph=graph, graphRandom=graphRandom, 
                                                      method=method,
                                                      ...,
                                                      FUN1=FUN, measure=measure,
                                                      verbose=verbose,
                                          rewire.w.type=rewire.w.type, 
                                           dist=dist,
                                          BPPARAM=BPPARAM)
    } else {
        
        if(type=="dependent")
        {
            print("Unweighted Network No Parallel Function")
            # No Parallel
            output <- robinRobustNoParallel(graph=graph, graphRandom= graphRandom, 
                                            method=method,
                                            ...,
                                            FUN=FUN, measure=measure,
                                            type=type, verbose=verbose) 
        }else{
            print("Unweighted Network Parallel Function")
            # Parallel version: 
            output <- robinRobustFast(graph=graph, graphRandom=graphRandom, 
                                      method=method,
                                      ...,
                                      FUN1=FUN, measure=measure,
                                      verbose=verbose,BPPARAM=BPPARAM)
        }
    }
    outputRobin <- c(output, model=methods, list(graph=graph))

class(outputRobin) <- "robin"
return(outputRobin)
}

####### GRAPH RANDOM #########
#' random
#'
#' @description This function randomly rewires the edges while preserving the original graph's 
#' degree distribution.
#' @param graph The output of prepGraph.
#' @param rewire.w.type for weighted graph. Option to rewire one of "Rewire",
#' "Shuffle","Garlaschelli","Sum" 
#' @param dist for weighted graph with "Garlaschelli" @rewire.w.type method. 
#' Option to rewire in a manner that retains overall graph weight regardless of 
#' distribution of edge weights. This option is invoked by putting any text into
#'  this field. Defaults to "Other". See \code{\link[perturbR]{rewireR}} for 
#'  details.
#' @param verbose flag for verbose output (default as FALSE)
#' 
#' @return An igraph object, a randomly rewired graph.
#' @import igraph
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)

 random <- function(graph, dist="Other", rewire.w.type="Rewire", verbose=FALSE)
{
    # Weigthed version
    if ( is_weighted(graph) )
    {
        graphRandom <- randomWeight(graph=graph,rewire.w.type=rewire.w.type,
                                    dist=dist, 
                                    verbose=verbose)
    }else{
        graphRandom <- randomNoW(graph=graph, verbose=verbose)
    }
    return(graphRandom)
}

