#' robinCompareFast
#'
#' @description This function compares two community detection algorithms.
#' Is the parallelized and faster version of \code{\link{robinCompare}}
#' @param graph The output of prepGraph.
#' @param method1 The first clustering method, one of "walktrap", 
#' "edgeBetweenness", "fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal".
#' @param method2 The second custering method one of "walktrap",
#' "edgeBetweenness","fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal".
#' @param FUN1 personal designed function when method1 is "others". 
#' see \code{\link{methodCommunity}}.
#' @param FUN2 personal designed function when method2 is "others". 
#' see \code{\link{methodCommunity}}.
#' @param measure The stability measure, one of "vi", "nmi", "split.join", 
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param weights This argument is not settable for "infomap" method.
#' @param steps This argument is settable only for "leadingEigen"and"walktrap" 
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method.
#' @param verbose flag for verbose output (default as TRUE).
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method 
#' - the matrix "Mean2" with the means of the procedure for the second method
#' 
#' @import igraph parallel
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' robinCompareFast(graph=graph, method1="louvain", 
#' method2="fastGreedy", measure="vi")

robinCompareFast <- function(graph, 
                         method1=c("walktrap", "edgeBetweenness", "fastGreedy",
                                   "leadingEigen","louvain","spinglass",
                                   "labelProp","infomap","optimal", "other"),
                         method2=c("walktrap", "edgeBetweenness", "fastGreedy",
                                   "leadingEigen","louvain","spinglass",
                                   "labelProp","infomap","optimal", "other"),
                         FUN1=NULL, FUN2=NULL,
                         measure= c("vi", "nmi","split.join", "adjusted.rand"),
                         directed=FALSE, weights=NULL, steps=4, 
                         spins=25, e.weights=NULL, v.weights=NULL, 
                         nb.trials=10, verbose=TRUE)
{   
    method1 <- match.arg(method1)
    method2 <- match.arg(method2)
    measure <-match.arg(measure)
    comReal1 <- membershipCommunities(graph=graph, method=method1,
                                      FUN=FUN1,
                                      directed=directed,
                                      weights=weights,
                                      steps=steps, 
                                      spins=spins, 
                                      e.weights=e.weights, 
                                      v.weights=v.weights, 
                                      nb.trials=nb.trials) 
    comReal2 <- membershipCommunities(graph=graph, method=method2,
                                      FUN=FUN2,
                                      directed=directed,
                                      weights=weights,
                                      steps=steps, 
                                      spins=spins, 
                                      e.weights=e.weights, 
                                      v.weights=v.weights, 
                                      nb.trials=nb.trials)
    de <- igraph::gsize(graph)
    N <- igraph::vcount(graph)
    Measure <- NULL
    vector1 <- NULL
    vector2 <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detecting robin method independent type, wait it can take
                    time it depends on the size of the network")
    vet1 <- seq(5, 60, 5) 
    vet <- round(vet1*de/100, 0)
    ncores <- parallel::detectCores(logical = FALSE) - 1
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl,varlist =c("graph","method1","method2","directed",
                                "weights","steps","spins", "e.weights", 
                                "v.weights", "nb.trials","measure","comReal1",
                                "comReal2","N","verbose","FUN1","FUN2"), 
                  envir=environment())
    zlist <- parallel::clusterApply(cl,vet, function(z) 
    {
        
        
        
        MeansList <- lapply(1:10, function(s)
        {
            
            graphRList <- igraph::rewire(graph, 
                                         with=igraph::keeping_degseq(loops=FALSE,
                                                                     niter=z))
            
            comr1 <- robin::membershipCommunities(graph=graphRList,
                                                  method=method1,
                                                  FUN=FUN1,
                                                  directed=directed,
                                                  weights=weights,
                                                  steps=steps, 
                                                  spins=spins, 
                                                  e.weights=e.weights, 
                                                  v.weights=v.weights, 
                                                  nb.trials=nb.trials)
            
            if(measure=="vi")
            {
                measure1 <- igraph::compare(comr1, comReal1, 
                                            method=measure)/log2(N)
            } else if(measure=="split.join")
            {
                measure1 <- igraph::compare(comr1, comReal1, 
                                            method=measure)/(2*N)
            }else{
                measure1 <- 1-(igraph::compare(comr1, comReal1, 
                                               method=measure))
                
            }
            
            
            
            comr2 <- robin::membershipCommunities(graph=graphRList, 
                                                  FUN=FUN2,
                                                  method=method2,
                                                  directed=directed,
                                                  weights=weights,
                                                  steps=steps, 
                                                  spins=spins, 
                                                  e.weights=e.weights, 
                                                  v.weights=v.weights, 
                                                  nb.trials=nb.trials)
            if(measure=="vi")
            {
                
                measure2 <-  igraph::compare(comr2, comReal2, 
                                             method=measure)/log2(N)
            } else if (measure=="split.join"){
                measure2 <- igraph::compare(comr2, comReal2, 
                                            method=measure)/(2*N)
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
    })
    stopCluster(cl)
    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))
    
    Measure1 <- cbind(rep(0, 10), Measure1)
    Measure2 <- cbind(rep(0, 10), Measure2)
    
    colnames(Measure1) <- nRewire 
    colnames(Measure2) <- nRewire 
    return(list(Mean1=Measure1,
                Mean2=Measure2))
}


