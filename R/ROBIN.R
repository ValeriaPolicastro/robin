###### PREPARATION GRAPH ########### 

#' prepGraph
#' 
#' @description This function reads graphs from a file and 
#' prepares them for the analysis.
#'
#' @param file The input file containing the graph.
#' @param file.format Character constant giving the file format. Edgelist, 
#' pajek, graphml, gml, ncol, lgl, dimacs, graphdb and igraph are
#' supported.
#' @param numbers A logical value indicating if the names of the nodes are 
#' values.This argument is settable for the edgelist format. 
#' The default is FALSE.
#' @param directed A logical value indicating if is a directed graph. The 
#' default is FALSE.
#' @param header A logical value indicating whether the file contains 
#' the names of the variables as its first line.This argument is settable 
#' @param verbose flag for verbose output (default as FALSE).
#' for the edgelist format.The default is FALSE.
#' @return An igraph object, which do not contain loop and multiple edges.
#' @import igraph
#' @importFrom utils read.table
#' @export
#'
#' @examples
#' #install.packages("robin")
#' 
#' #If there are problems with the installation try:
#' # if (!requireNamespace("BiocManager", quietly = TRUE))
#' #     install.packages("BiocManager")
#' # BiocManager::install("gprege")
#' # install.packages("robin")   
#'                      
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
prepGraph <- function(file, file.format=c("edgelist", "pajek", "ncol", "lgl",
                        "graphml", "dimacs", "graphdb", "gml", "dl","igraph"),
                        numbers=FALSE, directed=FALSE, header=FALSE, 
                        verbose=FALSE)
{ 
    file.format <- match.arg(file.format)
    if(verbose) cat("Detected file format: ", file.format, "\n")
    if (file.format =="igraph")
    { 
        net <- file
        ind <- igraph::V(net)[igraph::degree(net) == 0] #isolate node
        graph <- igraph::delete.vertices(net, ind)
        graph <- igraph::simplify(graph)
        
    }else if (file.format == "gml") {
        net <- igraph::read_graph(file=file, format=file.format)
        ind <- igraph::V(net)[degree(net) == 0] #isolate node
        graph <- igraph::delete.vertices(net, ind)
        graph <- igraph::simplify(graph)
    }else if((file.format == "edgelist") & (numbers == TRUE))
    {
        edge <- utils::read.table(file,colClasses = "character", quote="\"",
                                  header=header)
        edge <- as.matrix(edge)
        net <- igraph::graph_from_edgelist(edge, directed=directed)
        ind <- igraph::V(net)[degree(net) == 0] #isolate node
        graph <- igraph::delete.vertices(net, ind)
        graph <- igraph::simplify(graph)
    }else
    {
        net <- igraph::read_graph(file=file, format=file.format,
                                  directed=directed)
        ind <- igraph::V(net)[degree(net) == 0] #isolate node
        graph <- igraph::delete.vertices(net, ind)
        graph <- igraph::simplify(graph)
    }
    return(graph)
}


####### GRAPH RANDOM NO W #########
#' randomNoW
#'
#' @description This function randomly rewires the edges while preserving the original graph's 
#' degree distribution.
#' @param graph The output of prepGraph.
#' @param verbose flag for verbose output (default as FALSE)
#' 
#' @return An igraph object, a randomly rewired graph.
#' @import igraph
#' @keywords internal

randomNoW <- function(graph, verbose=FALSE)
{
    if(verbose) cat("Randomizing the graph edges.\n")
    z <- igraph::gsize(graph) ## number of edges
    graphRandom <- igraph::rewire(graph, 
                     with=igraph::keeping_degseq(loops=FALSE, niter=z))
    #rewiring for z all the edges
    return(graphRandom)
}


###### COMMUNITY METHOD ######    
#' methodCommunity
#' 
#' @description This function detects the community structure of a graph.
#' To detect the community structure the user can choose one of the methods implemented 
#' in igraph.
#' @param graph The output of prepGraph.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal", "leiden","other".
#' @param ... additional parameters to use with any of the previous described 
#' methods (see igraph package community detection methods for more details 
#' i.e. \link[igraph]{cluster_walktrap})
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param verbose flag for verbose output (default as FALSE)
#'
#' @return A Communities object.
#' @import igraph
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' methodCommunity (graph=graph, method="louvain") 
methodCommunity <- function(graph, 
                            method=c("walktrap", "edgeBetweenness", 
                                    "fastGreedy", "louvain", "spinglass", 
                                    "leadingEigen", "labelProp", "infomap",
                                    "optimal", "leiden", "other"), 
                                    ..., FUN=NULL, verbose=FALSE)
{   
    
    method <- match.arg(method)
    if(verbose) message("Applying community method ", method, "\n")
    dots <- list(...)
    communities <- switch(method, 
            optimal=igraph::cluster_optimal(graph, ...),
            louvain=igraph::cluster_louvain(graph=graph, ...),
            walktrap=igraph::cluster_walktrap(graph=graph, ...), 
            spinglass=igraph::cluster_spinglass(graph=graph, ...), 
            leadingEigen=igraph::cluster_leading_eigen(graph=graph, ...), 
            edgeBetweenness=igraph::cluster_edge_betweenness(graph=graph, ...), 
            fastGreedy=igraph::cluster_fast_greedy(graph=graph, ...), 
            labelProp=igraph::cluster_label_prop(graph=graph, ...), 
            infomap=igraph::cluster_infomap(graph=graph, ...),
            leiden=igraph::cluster_leiden(graph=graph, ...),
            other=FUN(graph, ...)
    )
    return(communities)
}

##### MEMBERSHIP COMMUNITIES ######    
#' membershipCommunities
#' 
#' @description This function computes the membership vector of the community 
#' structure. To detect the community structure the user can choose one of the methods implemented 
#' in igraph.
#' @param graph The output of prepGraph.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal", "leiden","other".
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param ... additional parameters to use with any of the previous described 
#' methods (see igraph package community detection methods for more details 
#' i.e. \link[igraph]{cluster_walktrap})
#'
#' @return Returns a numeric vector, one number for each vertex in the graph; 
#' the membership vector of the community structure.
#' @import igraph
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' membershipCommunities (graph=graph, method="louvain")

membershipCommunities <- function(graph,
                                 method=c("walktrap", "edgeBetweenness", 
                                        "fastGreedy", "louvain", "spinglass", 
                                        "leadingEigen", "labelProp", "infomap",
                                        "optimal", "leiden","other"), ...,
                                 FUN=NULL)
{
    method <- match.arg(method)
    members <- membership(methodCommunity(graph=graph, method=method, ...,
                                            FUN=FUN))
    
    return(members)
}


######### REWIRE COMPLETE ########
#' rewireCompl
#' 
#' @description rewires the graph, creates the communities and 
#' compares the communities through different measures.
#' 
#' @param data The output of prepGraph
#' @param number Number of rewiring trials to perform.
#' @param community Community to compare with.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap".
#' @param ... additional parameters to use with any of the previous described 
#' methods (see igraph package community detection methods for more details 
#' i.e. \link[igraph]{cluster_walktrap})
#' @param FUN see \code{\link{methodCommunity}}.
#' @param measure The measure for the comparison of the communities "vi", "nmi",
#' "split.join", "adjusted.rand"
#' @keywords internal
rewireCompl <- function(data, number, community, 
                        method=c("walktrap", "edgeBetweenness", 
                                 "fastGreedy", "louvain", "spinglass", 
                                 "leadingEigen", "labelProp", "infomap",
                                 "optimal", "leiden","other"), 
                        ..., 
                        measure=c("vi", "nmi","split.join", "adjusted.rand"),
                        FUN=NULL)
{
    method <- match.arg(method)
    measure <- match.arg(measure)
    # dots <- list(...)
    graphRewire <- igraph::rewire(data,
                            with=keeping_degseq(loops=FALSE, niter=number))
    comR <- membershipCommunities(graph=graphRewire, method=method, ..., FUN=FUN)
    Measure <- igraph::compare(community, comR, method=measure)
    output <- list(Measure=Measure, graphRewire=graphRewire)
    
   return(output)
}

######### REWIRE ONLY ###########
#' rewireOnl
#' @description makes the rewire function of igraph
#' @param data The output of prepGraph
#' @param number Number of rewiring trials to perform.
#' @keywords internal
rewireOnl <- function(data, number)
{
    graphRewire <- igraph::rewire(data, with=keeping_degseq(loops=FALSE,
                                              niter=number))
    return(graphRewire)
}


########  ROBIN PROCEDURE ROBUSTNESS#######
#' robinRobustNoParallel
#'
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
#' @param ... other parameter
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
#' robinRobustNoParallel(graph=graph, graphRandom=graphRandom, method="louvain",
#'     resolution=0.8, measure="vi", type="independent")
robinRobustNoParallel <- function(graph, graphRandom, 
                method=c("walktrap", "edgeBetweenness", 
                         "fastGreedy", "louvain", "spinglass", 
                         "leadingEigen", "labelProp", "infomap",
                         "optimal", "leiden", "other"),
                ...,
                FUN=NULL, measure= c("vi", "nmi","split.join", "adjusted.rand"),
                type=c("independent","dependent"),
                # directed=FALSE, weights=NULL, 
                # steps=4, spins=25, e.weights=NULL, v.weights=NULL, 
                # nb.trials=10, resolution=1, n_iterations=2,
                # objective_function = c("CPM", "modularity"), 
                verbose=TRUE)
{   
    measure <- match.arg(measure)
    type<- match.arg(type)
    method <- match.arg(method)
    # dots <- list(...)
     nrep <- 10
    comReal <- membershipCommunities(graph=graph, method=method, 
                                    ...=...,
                                    FUN=FUN) # real network
                                    # directed=directed,
                                    # weights=weights, 
                                    # steps=steps, 
                                    # spins=spins, 
                                    # e.weights=e.weights, 
                                    # v.weights=v.weights, 
                                    # nb.trials=nb.trials,
                                    # resolution=resolution,
                                    # objective_function = objective_function,
                                    # n_iterations=n_iterations
                                    # ) # real network

    comRandom <- membershipCommunities(graph=graphRandom, method=method,
                                       ...=..., FUN=FUN)
                                       #  directed=directed,
                                       #  weights=weights,
                                       #  steps=steps, 
                                       #  spins=spins, 
                                       #  e.weights=e.weights, 
                                       #  v.weights=v.weights, 
                                       #  nb.trials=nb.trials,
                                       # resolution=resolution,
                                       # objective_function = objective_function,
                                       # n_iterations=n_iterations) # random network
    #stopifnot(length(table(comRandom))>1)
    #if(length(table(comRandom))==1) {stop("Not random graph")}
    de <- igraph::gsize(graph)
    N <- igraph::vcount(graph)
    Measure <- NULL
    vector <- NULL
    vectRandom <- NULL
    graphRewireRandom <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0, 60, 5)
    if(verbose) cat("Detected robin method ", type, " type\n")
    #INDEPENDENT    
    if(type == "independent") 
    {
        #OUTPUT MATRIX
        measureReal <- matrix(0, nrep^2, length(nRewire))
        measureRandom <- matrix(0, nrep^2, length(nRewire))
        MeanRandom <- matrix(0, nrep, length(nRewire))
        Mean <- matrix(0, nrep, length(nRewire))
        vet1 <- seq(5, 60, 5) #each step 
        vet <- round(vet1*de/100, 0) #the numbers of edges to rewire
        #arrotonda a 0 cifre decimali
 
        for(z in vet)
        {
            count2 <- 0
            count <- count+1
            for(s in c(1:nrep))
            {
                count2 <- count2+1
                k <- 1
                #REAL
                Real <- rewireCompl(data=graph, 
                                    number=z, 
                                    community=comReal, 
                                    method=method,
                                    ...=..., FUN=FUN)
                                    # measure=measure,
                                    # directed=directed,
                                    # weights=weights,
                                    # steps=steps, 
                                    # spins=spins, 
                                    # e.weights=e.weights, 
                                    # v.weights=v.weights, 
                                    # nb.trials=nb.trials,
                                    # resolution=resolution,
                                    # objective_function = objective_function,
                                    # n_iterations=n_iterations)
                if (measure=="vi")
                {
                    vector[k] <- (Real$Measure)/log2(N)
                } else if(measure=="split.join"){
                     vector[k] <- (Real$Measure)/(2*N)
                } else {
                    vector[k] <- 1-(Real$Measure)
                }
                measureReal[count2, count] <- vector[k]
                graphRewire <- Real$graphRewire
                
                #RANDOM
                Random <- rewireCompl(data=graphRandom,
                                        number=z, 
                                        community=comRandom, 
                                        method=method,
                                        measure=measure, ...=..., FUN=FUN)
                                      #   directed=directed,
                                      #   weights=weights,
                                      #   steps=steps, 
                                      #   spins=spins, 
                                      #   e.weights=e.weights, 
                                      #   v.weights=v.weights, 
                                      #   nb.trials=nb.trials,
                                      # resolution=resolution,
                                      # objective_function = objective_function,
                                      # n_iterations=n_iterations)
                if (measure=="vi")
                {
                    vectRandom[k] <- (Random$Measure)/log2(N)
                } else if(measure=="split.join")
                {
                    vectRandom[k] <- (Random$Measure)/(2*N)
                } else {
                    vectRandom[k] <- 1-(Random$Measure)
                }
                measureRandom[count2, count] <- vectRandom[k]
                graphRewireRandom <- Random$graphRewire
                for(k in c(2:nrep))
                {
                    count2 <- count2+1
                    Real <- rewireCompl(data=graphRewire, 
                                        number=round(0.01*de), 
                                        community=comReal, 
                                        method=method,
                                        measure=measure, ...=..., FUN=FUN)
                                        # directed=directed,
                                        # weights=weights,
                                        # steps=steps, 
                                        # spins=spins, 
                                        # e.weights=e.weights, 
                                        # v.weights=v.weights, 
                                        # nb.trials=nb.trials,
                                        # resolution=resolution,
                                        # objective_function = objective_function,
                                        # n_iterations=n_iterations)
                    if (measure=="vi")
                    {
                        vector[k] <- (Real$Measure)/log2(N)
                    } else if(measure=="split.join")
                    {
                        vector[k] <- (Real$Measure)/(2*N)
                    } else {
                        vector[k] <- 1-(Real$Measure)
                    }
                    measureReal[count2, count] <- vector[k]
                    Random <- rewireCompl(data=graphRewireRandom, 
                                          number=round(0.01*de),
                                          community=comRandom,
                                          method=method,
                                          measure=measure, ...=..., FUN=FUN)
                                          # directed=directed,
                                          # weights=weights,
                                          # steps=steps, 
                                          # spins=spins, 
                                          # e.weights=e.weights, 
                                          # v.weights=v.weights, 
                                          # nb.trials=nb.trials,
                                          # resolution=resolution,
                                          # objective_function = objective_function,
                                          # n_iterations=n_iterations)
                    if (measure=="vi")
                    {
                        vectRandom[k] <- (Random$Measure)/log2(N)
                    } else if(measure=="split.join")
                    {
                        vectRandom[k] <- (Random$Measure)/(2*N)
                    } else {
                        vectRandom[k] <- 1-(Random$Measure)
                    }
                    measureRandom[count2, count] <- vectRandom[k]
                }
                MeanRandom[s, count] <- mean(vectRandom)
                Mean[s, count] <- mean(vector)
            }
            if(verbose) cat("Perturbed ", z, " edges\n")
        }
  #DEPENDENT 
    }else{
        z <- round((5*de)/100, 0) #the 5% of the edges
        measureReal <- rep(0, nrep^2)
        measureReal1 <- NULL
        measureRandom <- rep(0, nrep^2)
        measureRandom1 <- NULL
        MeanRandom <- rep(0, nrep)
        Mean <- rep(0, nrep)
        MeanRandom1 <- NULL
        Mean1 <- NULL
        diff <- NULL
        diffR <- NULL
        vet <- rep(z,(length(nRewire)-1))
        for(z in vet)
        {   
            count2 <- 0
            count <- count+1
            
            for(s in c(1:nrep))
                {
                count2 <- count2+1
                k <- 1
                ###REAL
                graphRewire <- rewireOnl(data=graph, number=z)
                graphRewire <-igraph::union(graphRewire, diff)
                comr <- membershipCommunities(graph=graphRewire,
                                        method=method, ...=..., FUN=FUN)
                                        # directed=directed,
                                        # weights=weights,
                                        # steps=steps, 
                                        # spins=spins, 
                                        # e.weights=e.weights, 
                                        # v.weights=v.weights, 
                                        # nb.trials=nb.trials,
                                        # FUN=FUN,
                                        # resolution=resolution,
                                        # objective_function = objective_function,
                                        # n_iterations=n_iterations)
                Measure <- igraph::compare(comReal, comr, method=measure)
                if (measure=="vi")
                {
                    vector[k] <- Measure/log2(N)
                } else if(measure=="split.join")
                {
                    vector[k] <- Measure/(2*N)
                } else {
                    vector[k] <- 1-Measure
                }
                measureReal1[count2] <- vector[k]
                diff <- igraph::difference(graph, graphRewire)
                
                ###RANDOM
                graphRewireRandom <- rewireOnl(data=graphRandom, number=z)
                graphRewireRandom <- igraph::union(graphRewireRandom, diffR)
                comr <- membershipCommunities(graph=graphRewireRandom,
                                        method=method, ...=..., FUN=FUN)
                                        # directed=directed,
                                        # weights=weights,
                                        # steps=steps, 
                                        # spins=spins, 
                                        # e.weights=e.weights, 
                                        # v.weights=v.weights, 
                                        # nb.trials=nb.trials,
                                        # FUN=FUN,
                                        # resolution=resolution,
                                        # objective_function = objective_function,
                                        # n_iterations=n_iterations)
                Measure <- igraph::compare(comRandom, comr,
                                           method=measure)
                if(measure=="vi")
                {
                    vectRandom[k] <- Measure/log2(N)
                } else if(measure=="split.join")
                {
                    vectRandom[k] <- Measure/(2*N)
                } else{
                    vectRandom[k] <- 1-Measure
                }
                
                measureRandom1[count2] <- vectRandom[k]
                diffR <- igraph::difference(graphRandom, graphRewireRandom)
                
                for(k in c(2:nrep)) 
                {
                    count2 <- count2+1
                    ##REAL
                    Real <- rewireCompl(data=graphRewire, number=round(0.01*de),
                                        method=method,
                                        measure=measure, ...=..., FUN=FUN)
                                        # community=comReal,
                                        # directed=directed,
                                        # weights=weights,
                                        # steps=steps, 
                                        # spins=spins, 
                                        # e.weights=e.weights, 
                                        # v.weights=v.weights, 
                                        # nb.trials=nb.trials,
                                        # resolution=resolution,
                                        # objective_function = objective_function,
                                        # n_iterations=n_iterations)
                    if(measure=="vi")
                    {
                        vector[k] <- (Real$Measure)/log2(N)
                    } else if(measure=="split.join")
                    {
                        vector[k] <- (Real$Measure)/(2*N)
                    } else{
                        vector[k] <- 1-(Real$Measure)
                    }
                    measureReal1[count2] <- vector[k]
                    
                    ## RANDOM
                    Random <- rewireCompl(data=graphRewireRandom,
                                            method=method,
                                            measure=measure,
                                            number=round(0.01*de), 
                                            ...=..., FUN=FUN)
                                          #   community=comRandom,
                                          #   directed=directed,
                                          #   weights=weights,
                                          #   steps=steps, 
                                          #   spins=spins, 
                                          #   e.weights=e.weights, 
                                          #   v.weights=v.weights, 
                                          #   nb.trials=nb.trials,
                                          # resolution=resolution,
                                          # objective_function = objective_function,
                                          # n_iterations=n_iterations)
                    if(measure=="vi")
                    {
                        vectRandom[k] <- (Random$Measure)/log2(N)
                    } else if(measure=="split.join")
                    {
                        vectRandom[k] <- (Random$Measure)/(2*N)
                    } else{
                        vectRandom[k] <- 1-(Random$Measure)
                    }
                    measureRandom1[count2] <- vectRandom[k]
                }
                Mean1[s] <- mean(measureReal1)
                MeanRandom1[s] <- mean(measureRandom1)  
            }
            graph <- igraph::intersection(graph, graphRewire)
            graphRandom <- igraph::intersection(graphRandom, graphRewireRandom)
            measureRandom <- cbind(measureRandom, measureRandom1)
            measureReal <- cbind(measureReal, measureReal1)
            Mean <- cbind(Mean, Mean1)
            MeanRandom <- cbind(MeanRandom, MeanRandom1)
            z1 <- igraph::gsize(graph)
            #print(z1)
            if(verbose) cat("Perturbed ", z, " edges\n")
            }
    }
    colnames(measureRandom) <- nRewire
    colnames(measureReal) <- nRewire
    #the matrices "measureReal" and "measureRandom" with the 
    #measures calculated at each step of the procedure, respectively for the real and 
    #the random graph.
    colnames(MeanRandom) <- nRewire
    colnames(Mean) <- nRewire
    output <- list( Mean=Mean,
                    MeanRandom=MeanRandom
                    )
      return(output)

}


############### COMPARISON DIFFERENT METHODS ##########
#' robinCompareNoParallel
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
#' @param verbose flag for verbose output (default as TRUE).
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method 
#' - the matrix "Mean2" with the means of the procedure for the second method
#' 
#' @import igraph
#' @export 
#' @keywords internal
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' robinCompareNoParallel(graph=graph, method1="louvain", args1 = list(resolution=0.8),
#'             method2="leiden", args2=list(objective_function ="modularity"), 
#'             measure="vi", type="independent")
robinCompareNoParallel <- function(graph, 
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
                         type=c("independent", "dependent"), verbose=TRUE)
{
    method1 <- match.arg(method1)
    method2 <- match.arg(method2)
    type <- match.arg(type)
    measure <- match.arg(measure)
    nrep <- 10
    N <- igraph::vcount(graph)
    args11 <- c(list(graph=graph), method=method1, FUN=FUN1, args1)
    args21 <- c(list(graph=graph), method=method2, FUN=FUN2, args2)
    comReal1 <- do.call(membershipCommunities, args11)
    comReal2 <- do.call(membershipCommunities, args21)
    
    de <- igraph::gsize(graph)
    Measure <- NULL
    vector1 <- NULL
    vector2 <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detected robin method ", type, " type\n")
    #independent
    if(type == "independent") 
    {
        measureReal1 <- matrix(0, nrep^2, length(nRewire))
        measureReal2 <- matrix(0, nrep^2, length(nRewire))
        Mean1 <- matrix(0, nrep, length(nRewire))
        Mean2 <- matrix(0, nrep, length(nRewire))
        vet1 <- seq(5, 60, 5) 
        vet <- round(vet1*de/100, 0)
        
        for(z in vet)
        {
            count2 <- 0
            count <- count+1
            for(s in c(1:nrep))
            {
                count2 <- count2+1
                k <- 1
                graphRewire <- rewireOnl(data=graph, number=z)
                args11 <- c(list(graph=graphRewire), method=method1, FUN=FUN1, args1)
                args21 <- c(list(graph=graphRewire), method=method2, FUN=FUN2, args2)
                comr1 <- do.call(membershipCommunities, args11)
                comr2 <- do.call(membershipCommunities, args21)
                if (measure=="vi")
                {
                    vector1[k] <- igraph::compare(comr1, comReal1, 
                                                  method=measure)/log2(N)
                    vector2[k] <- igraph::compare(comr2, comReal2, 
                                                  method=measure)/log2(N)
                } else if (measure=="split.join")
                {
                    vector1[k] <- igraph::compare(comr1, comReal1, 
                                                  method=measure)/(2*N)
                    vector2[k] <- igraph::compare(comr2, comReal2, 
                                                  method=measure)/(2*N)
                } else {
                    vector1[k] <- 1-(igraph::compare(comr1, comReal1, 
                                                     method=measure))
                    vector2[k] <- 1-(igraph::compare(comr2, comReal2, 
                                                     method=measure))
                }
                measureReal1[count2, count] <- vector1[k]
                measureReal2[count2, count] <- vector2[k]
                
                for(k in c(2:nrep))
                {
                    count2 <- count2+1
                    graphRewire <- rewireOnl(data=graphRewire,
                                             number=round(0.01*z))
                    args11 <- c(list(graph=graphRewire), method=method1, FUN=FUN1, args1)
                    args21 <- c(list(graph=graphRewire), method=method2, FUN=FUN2, args2)
                    comr1 <- do.call(membershipCommunities, args11)
                    comr2 <- do.call(membershipCommunities, args21)
                   
                    if(measure=="vi")
                    {
                        vector1[k] <- igraph::compare(comr1, comReal1, 
                                                      method=measure)/log2(N)
                        vector2[k] <- igraph::compare(comr2, comReal2, 
                                                      method=measure)/log2(N)
                    } else if(measure=="split.join")
                    {
                        vector1[k] <- igraph::compare(comr1, comReal1, 
                                                      method=measure)/(2*N)
                        vector2[k] <- igraph::compare(comr2, comReal2, 
                                                      method=measure)/(2*N)
                    } else{
                        vector1[k] <- 1-(igraph::compare(comr1, comReal1, 
                                                         method=measure))
                        vector2[k] <- 1-(igraph::compare(comr2, comReal2, 
                                                         method=measure))
                    }
                    
                    measureReal1[count2, count] <- vector1[k]
                    measureReal2[count2, count] <- vector2[k]
                    
                    
                }
                Mean1[s, count] <- mean(vector1)
                Mean2[s, count] <- mean(vector2)
            }
            if(verbose) cat("Perturbed ", z, " edges\n")
        }
        #dependent    
    }else{
        z <- round((5*de)/100, 0)
        measureReal1 <- rep(0, nrep^2)
        measureReal11 <- NULL
        R11<-NULL
        measureReal2 <- rep(0, nrep^2)
        measureReal22 <- NULL
        R22 <- NULL
        Mean1 <- rep(0, nrep)
        Mean2 <-rep(0, nrep)
        Mean11 <- NULL
        Mean22 <- NULL
        diff <- NULL
        vet<-rep(z,(length(nRewire)-1))
        for(z in vet)
        {   
            count2 <- 0
            count <- count+1
            for(s in c(1:nrep))
            {
                count2 <- count2+1
                k <- 1
                graphRewire <- rewireOnl(data=graph, number=z)
                graphRewire <- igraph::union(graphRewire, diff)
                args11 <- c(list(graph=graphRewire), method=method1, FUN=FUN1, args1)
                args21 <- c(list(graph=graphRewire), method=method2, FUN=FUN2, args2)
                comr1 <- do.call(membershipCommunities, args11)
                comr2 <- do.call(membershipCommunities, args21)
                
                if(measure=="vi")
                {
                    vector1[k] <- igraph::compare(comr1, comReal1, 
                                                  method= measure)/log2(N)
                    vector2[k] <- igraph::compare(comr2, comReal2, 
                                                  method= measure)/log2(N)
                } else if(measure=="split.join")
                {
                    vector1[k] <- igraph::compare(comr1, comReal1, 
                                                  method= measure)/(2*N)
                    vector2[k] <- igraph::compare(comr2, comReal2, 
                                                  method= measure)/(2*N)
                } else{
                    vector1[k] <- 1-(igraph::compare(comr1, comReal1, 
                                                     method= measure))
                    vector2[k] <- 1-(igraph::compare(comr2, comReal2, 
                                                     method= measure))
                }
                
                measureReal11[count2] <- vector1[k]
                measureReal22[count2] <- vector2[k]
                diff <- igraph::difference(graph, graphRewire)
                
                for(k in c(2:nrep)) 
                {
                    count2 <- count2+1
                    graphRewire <- rewireOnl(data=graphRewire,
                                             number=round(0.01*z))
                    args11 <- c(list(graph=graphRewire), method=method1, FUN=FUN1, args1)
                    args21 <- c(list(graph=graphRewire), method=method2, FUN=FUN2, args2)
                    comr1 <- do.call(membershipCommunities, args11)
                    comr2 <- do.call(membershipCommunities, args21)
                    
                    if(measure=="vi")
                    {
                        vector1[k] <- igraph::compare(comr1, comReal1, 
                                                      method=measure)/log2(N)
                        vector2[k] <- igraph::compare(comr2, comReal2, 
                                                      method=measure)/log2(N)
                    } else  if(measure=="split.join")
                    {
                        vector1[k] <- igraph::compare(comr1, comReal1, 
                                                      method=measure)/(2*N)
                        vector2[k] <- igraph::compare(comr2, comReal2, 
                                                      method=measure)/(2*N)
                    } else{
                        vector1[k] <- 1-(igraph::compare(comr1, comReal1, 
                                                         method=measure))
                        vector2[k] <- 1-(igraph::compare(comr2, comReal2, 
                                                         method=measure))
                    }
                    measureReal11[count2] <- vector1[k]
                    measureReal22[count2] <- vector2[k]
                    
                }
                Mean11[s] <- mean(measureReal11)
                Mean22[s] <- mean(measureReal22)
                
            }
            graph <- igraph::intersection(graph, graphRewire)
            Mean1 <-cbind(Mean1,Mean11)
            Mean2 <-cbind(Mean2,Mean22)
            z1 <- igraph::gsize(graph)
            if(verbose) cat("Perturbed ", z, " edges\n")
        }
    }
    colnames(Mean1) <- nRewire 
    colnames(Mean2) <- nRewire 
    output <- list(Mean1=Mean1,
                   Mean2=Mean2)
    return(output)
}



########################### TEST ROBIN ###################

################ CREATE ITPSpline ##################

#' @title createITPSplineResult
#' 
#' @description creates an fdatest::ITP2 class object 
#' 
#' @param graph The output of prepGraph.
#' @param model1 The Mean output of the robinRobust function (or the Mean1 
#' output of the comparison function).
#' @param model2 The MeanRandom output of the robinRobust function (or the 
#' Mean2 output of the comparison function).
#' @param muParam the mu parameter for ITP2bspline (default 0).
#' @param orderParam the order parameter for ITP2bspline (default 4).
#' @param nKnots the nknots parameter for ITP2bspline (default 7).
#' @param BParam the B parameter for ITP2bspline (default 10000).
#' @param isPaired the paired parameter for ITP2bspline (default TRUE).
#' 
#' @keywords internal
#' 
#' 

createITPSplineResult <- function(graph, model1, model2,
                                muParam=0, orderParam=4, nKnots=7, 
                                BParam=10000, isPaired=TRUE) 
{

    modeled1 <- as.matrix(model1)
    modeled2 <- as.matrix(model2)   
    
   
    
    ## Two populations Interval Testing Procedure with B-spline basis
    ITPresult <- fdatest::ITP2bspline(modeled1, modeled2, mu=muParam, 
                                    order=orderParam, nknots=nKnots, 
                                    B=BParam, paired=isPaired)
   
    return(ITPresult)
}

######################GAUSSIAN PROCESS#################

#' robinGPTest
#'
#' @description This function implements the GP testing procedure and calculates the 
#' Bayes factor.
#' @param x A robin class object. The output of the functions:  
#' \code{\link{robinRobust}} and \code{\link{robinCompare}}.
#' @param verbose flag for verbose output (default as FALSE).
#' 
#' @return A numeric value, the Bayes factor
#' @importFrom stats sd var
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' comp <- robinCompare(graph=graph, method1="fastGreedy",method2="infomap")
#' \donttest{robinGPTest(comp)}
robinGPTest <- function(x, verbose=FALSE)
{ 
   
    if (length(x$Mean1)==0)
    {
        model1 <- x$Mean
        model2 <- x$MeanRandom 
    }else{
        model1 <- x$Mean1
        model2 <- x$Mean2 
    }
    ratios <- log2((model1+0.001)/(model2+0.001))
   #rapporto tra la media delle misure tra il modello reale e quello perturbato 
   #e la media delle distanze tra il random e la sua perturbazione
   res <- as.vector(ratios)

   nRewire <- as.numeric(colnames(model1)) # generico anche se faccio un rewire fino al 40%
   nrep <- dim(model1)[1]
   names(res) <- rep(nRewire, each=nrep)

   ratio <- t(res)#la trasposta del rapporto
   gpregeOptions <- list(indexRange=(1:2), explore=FALSE,
                         exhaustPlotRes=30, exhaustPlotLevels=10,
                         exhaustPlotMaxWidth=100, iters=100,
                         labels=rep(FALSE,2), display=FALSE)

    if(verbose) cat("Computing Gaussian Process Testing.\n")
    MA <- as.matrix(ratio[,c(1:dim(ratio)[2])])
    MA <- t(MA)
    vt <- unique(colnames(MA))
    ntimes <- length(vt)
    stdv <- NULL
    varv <- NULL
    for (i in c(1:ntimes)){    #ntime number of percentuage of rewire
        ind <- which(colnames(MA)==vt[i])
        stdv[i] <- stats::sd(MA[1,ind])
        varv[i] <- stats::var(MA[1,ind])
    }
    GlobalVar <- var(MA[1,])
    SigNoise <- mean(varv)/GlobalVar
    if (SigNoise>1) SigNoise <- 1
    sigmaest <- 1-SigNoise
    ## [inverse-lengthscale percent-signal-variance percent-noise-variance]
    gpregeOptions$inithypers <- matrix( c(
        1/1000,	0,	1
        ,1/ntimes,	sigmaest, SigNoise 
    ), ncol=3, byrow=TRUE)
    dvet <- data.matrix(as.numeric(colnames(ratio)))
    dd <- t(data.matrix(as.numeric(ratio)))
    rownames(dd) <- 'Measure'
    colnames(dd) <- dvet
    datadum <- rbind(dd, dd)
    gpregeOutput <- .gprege(data=datadum, inputs=dvet,
                                  gpregeOptions=gpregeOptions)
    bf <- gpregeOutput$rankingScores[1]

    return(Bayes_Factor=bf)
  
 }


################ FDA TEST ##############################

#' robinFDATest
#'
#'@description The function implements the Interval Testing Procedure to 
#'test the difference between two curves.
#'
#' @param x A robin class object. The output of the functions:  
#' \code{\link{robinRobust}} and \code{\link{robinCompare}}.
#' @param verbose flag for verbose output (default as FALSE).
#' 
#' @return Two plots: the fitted curves and the adjusted p-values. A vector of the adjusted p-values. 
#' @import qpdf igraph ggplot2 fdatest graphics
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' comp <- robinCompare(graph=graph, method1="fastGreedy",method2="infomap")
#' \donttest{robinFDATest(comp)}
robinFDATest <- function(x, verbose=FALSE)
{
    if(verbose) cat("Computing Interval testing procedure.\n")
    
    graph <- x$graph
    legend <- c(x$model1,x$model2)
    if (length(x$Mean1)==0)
    {
        model1 <- x$Mean
        model2 <- x$MeanRandom 
    }else{
        model1 <- x$Mean1
        model2 <- x$Mean2 
    }
    
    
    
    object <- createITPSplineResult(graph, model1, model2)
    #Functional Data plot
    J <- dim(object$data.eval)[2]
    xmin <- 0
    xmax <- 0.6
    Abscissa <- seq(xmin,xmax,len=J)
    
    model1 <- cbind(as.numeric(as.vector(t(object$data.eval[1:10,]))))
    model2 <- cbind(as.numeric(as.vector(t(object$data.eval[11:20,]))))
    
    measures <- rbind(model1, model2)
    model <- c(rep(legend[1],each=10000),rep(legend[2],each=10000))
    percPert <- as.numeric(rep(Abscissa, times = 10))
    s <- c(rep(1:10,each=1000),rep(11:20,each=1000))
    dataFrame <- data.frame(measures,model,s,percPert)
    plot1 <- ggplot2::ggplot(dataFrame, ggplot2::aes(x=as.numeric(percPert),
                 y=as.numeric(measures), color= model, group=s)) +
             ggplot2::geom_line() +
             ggplot2::xlab("Percentage of perturbation") +
             ggplot2::ylab("Measure")+
             ggplot2::ggtitle("Functional Data Analysis")+
             ggplot2::scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6))
      
    
     
  
     #P value plot
     p <- length(object$pval)
     xmin <- 0
     xmax <- 0.6
     abscissa.pval <- rep(seq(xmin,xmax,len=p),time=2)
     pvalue <- c(object$pval,object$corrected.pval)
     type <- c(rep("pvalue",p),rep("pvalue.adj",p))
     PdataFrame <- data.frame(cbind(abscissa.pval,pvalue,type))
     plot2 <- ggplot2::ggplot(PdataFrame, ggplot2::aes(x=as.numeric(abscissa.pval),
                                                    y=as.numeric(pvalue), color= type)) +
              ggplot2::geom_point() +
              ggplot2::xlab("Percentage of perturbation") +
              ggplot2::ylab("p_value")+
              ggplot2::ggtitle("P-values")+
              ggplot2::scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6))+
              ggplot2::geom_hline(yintercept = 0.05,color = "red")+
              ggplot2::lims(y=c(0,1))
     
    
    plot <- gridExtra::grid.arrange(plot1,plot2, ncol=2)
    print(plot)
   
    adj.pvalue <- object$corrected.pval
    pvalue <- object$pval
    output <- list(adj.pvalue=adj.pvalue,
                 pvalues=pvalue)
    
    return(output)
}  

########### AREA UNDER THE CURVE    ##############
#' robinAUC
#'
#' @description This function calculates the area under two curves with a spline approach. 
#' @param x A robin class object. The output of the functions:  
#' \code{\link{robinRobust}} and \code{\link{robinCompare}}.
#' @param verbose flag for verbose output (default as FALSE).
#' 
#' @return A list
#' @importFrom DescTools AUC
#' @importFrom igraph vcount
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' proc <- robinRobust(graph=graph, graphRandom=graphRandom, method="louvain",
#' measure="vi")
#' robinAUC(proc)
robinAUC <- function( x, verbose=FALSE)
{
    if(verbose) cat("Computing area under the curve (AUC).\n")
    
    if (length(x$Mean1)==0)
    {
        model1 <- x$Mean
        model2 <- x$MeanRandom 
    }else{
        model1 <- x$Mean1
        model2 <- x$Mean2 
    }
    
        mvimeanmodel1 <- cbind(as.vector((apply(model1, 2, mean))))
        mvimeanmodel2 <- cbind(as.vector((apply(model2, 2, mean))))
    
    area1 <- DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel1, 
                          method ="spline")
    area2 <- DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel2, 
                          method ="spline")
    
    output <- c(area1,area2)
    names(output) <- c(paste("Area",x$model1),paste("Area", x$model2))
    
return(output)
}

