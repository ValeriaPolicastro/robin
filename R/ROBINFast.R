#' robinCompareFast
#'
#' @description This function compares two community detection algorithms.
#' Is the parallelized and faster version of \code{\link{robinCompare}}
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
#' @param type The type of robin construction, dependent or independent.
#' @param ncores number of CPU cores to use.(default is 2) For a faster 
#' execution we suggest to use ncores=(parallel::detectCores(logical = FALSE)-1) 
#' @param verbose flag for verbose output (default as TRUE).
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method 
#' - the matrix "Mean2" with the means of the procedure for the second method
#' 
#' @import igraph parallel
#' @keywords internal
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' robinCompareFast(graph=graph, method1="louvain", args1 = list(resolution=0.8),
#'             method2="leiden", args2=list(objective_function ="modularity"),
#'             ncores=2)

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
                         ncores=2,
                         FUN1=NULL, FUN2=NULL,
                         verbose=TRUE)
{   
    method1 <- match.arg(method1)
    method2 <- match.arg(method2)
    measure <- match.arg(measure)
    args11 <- c(list(graph=graph), method=method1, FUN=FUN1, args1)
    args21 <- c(list(graph=graph), method=method2, FUN=FUN2, args2)
    comReal1 <- do.call(robin::membershipCommunities, args11)
    comReal2 <- do.call(robin::membershipCommunities, args21)
    
    # comReal1 <- membershipCommunities(graph=graph, method=method1,
    #                                   FUN=FUN1,
    #                                   directed=directed,
    #                                   weights=weights,
    #                                   steps=steps, 
    #                                   spins=spins, 
    #                                   e.weights=e.weights, 
    #                                   v.weights=v.weights, 
    #                                   nb.trials=nb.trials,
    #                                   resolution=resolution,
    #                                   objective_function = objective_function,
    #                                   n_iterations=n_iterations) 
    # comReal2 <- membershipCommunities(graph=graph, method=method2,
    #                                   FUN=FUN2,
    #                                   directed=directed,
    #                                   weights=weights,
    #                                   steps=steps, 
    #                                   spins=spins, 
    #                                   e.weights=e.weights, 
    #                                   v.weights=v.weights, 
    #                                   nb.trials=nb.trials,
    #                                   resolution=resolution,
    #                                   objective_function = objective_function,
    #                                   n_iterations=n_iterations)
    de <- igraph::gsize(graph)
    N <- igraph::vcount(graph)
   # Measure <- NULL
     # vector1 <- NULL
   # vector2 <- NULL
    # graphRewire <- NULL
   # count <- 1
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detecting robin method independent parallelized type, wait it can take time it depends on the size of the network.\n")
    vet1 <- seq(5, 60, 5) 
    vet <- round(vet1*de/100, 0)
    cl <- parallel::makeCluster(ncores)
    
    #varli <- c(list(graph=graph),method1=method1, method2=method2, FUN1=FUN1, 
    #             FUN2=FUN2)
    #varli <- c(names(args1),
    #           names(args2),
    #           names(varli))
    #args11 <- c(cl=cl, varlist=varlist, envir=environment())
    #do.call(parallel::clusterExport, args11)
    
     parallel::clusterExport(cl,varlist=c("graph", 
                                          "method1", 
                                         "method2", 
                                          "FUN1",
                                         "FUN2",
                                         "args1",
                                         "args2",
                                         "comReal1",
                                         "comReal2"), envir=environment())
    #print(names(varlist))
    zlist <- parallel::clusterApply(cl, vet, function(z) 
    {
        
        
        
        MeansList <- lapply(1:10, function(s)
        {
            
            graphRList <- igraph::rewire(graph, 
                                with=igraph::keeping_degseq(loops=FALSE,
                                                             niter=z))
            
            # comr1 <- robin::membershipCommunities(graph=graphRList,
            #                                       method=method1,
            #                                       FUN=FUN1,
            #                                       directed=directed,
            #                                       weights=weights,
            #                                       steps=steps, 
            #                                       spins=spins, 
            #                                       e.weights=e.weights, 
            #                                       v.weights=v.weights, 
            #                                       nb.trials=nb.trials,
            #                                       resolution=resolution,
            #                                       objective_function = objective_function,
            #                                       n_iterations=n_iterations)
           
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
            # comr2 <- robin::membershipCommunities(graph=graphRList, 
            #                                       FUN=FUN2,
            #                                       method=method2,
            #                                       directed=directed,
            #                                       weights=weights,
            #                                       steps=steps, 
            #                                       spins=spins, 
            #                                       e.weights=e.weights, 
            #                                       v.weights=v.weights, 
            #                                       nb.trials=nb.trials,
            #                                       resolution=resolution,
            #                                       objective_function = objective_function,
            #                                       n_iterations=n_iterations)
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
    })
    parallel::stopCluster(cl)
    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))
    
    Measure1 <- cbind(rep(0, 10), Measure1)
    Measure2 <- cbind(rep(0, 10), Measure2)
    
    colnames(Measure1) <- nRewire 
    colnames(Measure2) <- nRewire 
    return(list(Mean1=Measure1,
                Mean2=Measure2))
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
#' @param ... 
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param ... other parameter
#' @param verbose flag for verbose output (default as TRUE)
#' @param BPPARAM the BiocParallel object of class \code{bpparamClass} that 
#' specifies the back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean" with the means of the procedure for the graph
#' - the matrix "MeanRandom" with the means of the procedure for the random graph. 
#' 
#' @import igraph
#' @importFrom BiocParallel bplapply bpparam
#' @export
#' @keywords internal
#' @example 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' robinRobustFast(graph=graph, graphRandom=graphRandom, method="leiden",
#'     resolution_parameter = 1, measure="vi")
robinRobustFast <- function(graph, graphRandom, 
                            method=c("walktrap", "edgeBetweenness", 
                                     "fastGreedy", "louvain", "spinglass", 
                                     "leadingEigen", "labelProp", "infomap",
                                     "optimal", "leiden", "other"),
                            ...,
                            FUN=NULL, measure= c("vi", "nmi", "split.join", 
                                                 "adjusted.rand"),
                            verbose=TRUE, BPPARAM=BiocParallel::bpparam())
{
    method <- match.arg(method)
    measure <- match.arg(measure)
    comReal1 <- membershipCommunities(graph=graph, method=method,
                                      FUN=FUN, ...=...) 
    comReal2 <- membershipCommunities(graph=graphRandom, method=method,
                                      FUN=FUN, ...=...)
    de <- igraph::gsize(graph)
    N <- igraph::vcount(graph)
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detecting robin method independent type, wait it can take time it depends on the size of the network.\n")
    vet1 <- seq(5, 60, 5) 
    vet <- round(vet1*de/100, 0)
    ## read ... and unpack -> create a list of character strings and pass it 
    ## in varlist
    # varlist <- list(...)
    # varlist <- c(varlist, "graph","method","comReal1",
    #              "comReal2","N","verbose","FUN", 
    #              "measure")
    # cl <- parallel::makeCluster(ncores)
    # parallel::clusterExport(cl, varlist=c("...", "graph","method","comReal1",
    #                                      "comReal2","N","verbose","FUN", 
    #                                      "measure"),
    #                        envir=environment())
    # zlist <- parallel::clusterApply(cl, vet, function(z) 
    
    parfunct <- function(z, graph, method, comReal1, comReal2, N, 
                         measure, ...)
    {
        print(list(...))
        MeansList <- lapply(1:10, function(s)
        {
            
            graphRList <- igraph::rewire(graph, 
                                         with=igraph::keeping_degseq(loops=FALSE,
                                                                     niter=z))
            
            comr1 <- robin::membershipCommunities(graph=graphRList,
                                                  method=method,
                                                  FUN=FUN,
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
                                                  FUN=FUN,
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
    
    zlist <- BiocParallel::bplapply(vet, parfunct, graph=graph, measure=measure,
            method=method, comReal1=comReal1, N=N, comReal2=comReal2, ...=...,
            BPPARAM=BPPARAM)
    
    # parallel::stopCluster(cl)
    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))
    
    Measure1 <- cbind(rep(0, 10), Measure1)
    Measure2 <- cbind(rep(0, 10), Measure2)
    
    colnames(Measure1) <- nRewire 
    colnames(Measure2) <- nRewire 
    return(list(Mean=Measure1,
                MeanRandom=Measure2))
}


