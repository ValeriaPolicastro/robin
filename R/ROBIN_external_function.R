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
#' @param type The type of robin construction, parallel, dependent or 
#' independent (default parallel).
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
                          measure= c("vi", "nmi","split.join", "adjusted.rand"),
                          ncores=2,
                          type=c("parallel", "independent", "dependent"),
                          verbose=TRUE)
{
    type <- match.arg(type)
    if(type=="parallel"){
        
     output <- robinCompareFast(graph=graph, method1=method1, args1 = args1,
                         method2=method2, args2=args2, measure=measure, 
                         ncores=ncores)
        
    }
    
    else 
    {
        output <- robinCompareNoParallel(graph=graph, method1=method1, args1 = args1,
                            method2=method2, args2=args2, measure=measure, 
                           type=type) 
    }
   
    
   # class(output) <- "robin"
    return(output)
    
}
           