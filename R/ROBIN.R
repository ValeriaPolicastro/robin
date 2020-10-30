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
        graph <- igraph::simplify(file) 
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
        graph <- igraph::graph_from_edgelist(edge, directed=directed)
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


####### GRAPH RANDOM #########
#' random
#'
#' @description This function randomly rewires the edges while preserving the original graph's 
#' degree distribution.
#' @param graph The output of prepGraph.
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
random <- function(graph, verbose=FALSE)
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
#' "optimal", "other".
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param weights  Optional positive weight vector. If the graph has a weight 
#' edge attribute, then this is used by default. Supply NA here if the graph 
#' has a weight edge attribute, but you want to ignore it. Larger edge weights
#' correspond to stronger connections. This argument is not settable for 
#' "infomap" method.
#' @param steps The number of steps to take, this is actually the number of 
#' tries to make a step. It is not a particularly useful parameter. This 
#' argument is settable only for "leadingEigen" and "walktrap" method.
#' @param spins Integer constant, the number of spins to use. This is the upper 
#' limit for the number of communities. It is not a problem to supply a 
#' (reasonably) big number here, in which case some spin states will be 
#' unpopulated. This argument is settable only for "spinglass" method.
#' @param e.weights If not NULL, then a numeric vector of edge weights. 
#' The length must match the number of edges in the graph. By default the 
#' ‘weight’ edge attribute is used as weights. If it is not present, then all
#' edges are considered to have the same weight. Larger edge weights correspond 
#' to stronger connections. This argument is settable only for "infomap"
#'  method.
#' @param v.weights If not NULL, then a numeric vector of vertex weights. The
#' length must match the number of vertices in the graph. By default the 
#' ‘weight’ vertex attribute is used as weights. If it is not present, then all
#' vertices are considered to have the same weight. A larger vertex weight means
#' a larger probability that the random surfer jumps to that vertex. This 
#' argument is settable only for "infomap" method.
#' @param nb.trials The number of attempts to partition the network (can be any
#' integer value equal or larger than 1). This argument is settable only for
#' "infomap" method.
#' @param directed Logical constant, whether to calculate directed edge 
#' betweenness for directed graphs. This argument is settable only for 
#' "edgeBetweenness" method.
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
                                    "optimal", "other"),
                            FUN=NULL, directed=FALSE, weights=NULL, steps=4, 
                            spins=25, e.weights=NULL, v.weights=NULL, 
                            nb.trials=10, verbose=FALSE)
{   
    
    method <- match.arg(method)
    if(verbose) cat("Applying community method ", method, "\n")
    if(is.null(weights) &
       (sum(method %in% c("walktrap", "edgeBetweenness", "fastGreedy") == 1 )))
    {
        weights <- E(graph)$weight
    }
    
    if((steps == 4) & (method == "leadingEigen"))
    {
        steps  <- -1
    }
    communities <- switch(method, 
            optimal=igraph::cluster_optimal(graph, weights = weights),
            louvain=igraph::cluster_louvain(graph=graph, weights=weights),
            walktrap=igraph::cluster_walktrap(graph=graph, weights=weights, 
                                    steps=steps), 
            spinglass=igraph::cluster_spinglass(graph=graph, weights=weights, 
                                        spins=spins), 
            leadingEigen=igraph::cluster_leading_eigen(graph=graph, 
                                            steps=steps, weights=weights,
                                            options=list(maxiter=1000000)), 
            edgeBetweenness=igraph::cluster_edge_betweenness(graph=graph, 
                                                    weights=weights,
                                                    directed=directed), 
            fastGreedy=igraph::cluster_fast_greedy(graph=graph, 
                                                   weights=weights), 
            labelProp=igraph::cluster_label_prop(graph=graph, weights=weights), 
            infomap=igraph::cluster_infomap(graph=graph, e.weights=e.weights, 
                                v.weights=v.weights, nb.trials=nb.trials),
            other=FUN(graph, weights)
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
#' "optimal".
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param weights  Optional positive weight vector. If the graph has a weight 
#' edge attribute, then this is used by default. Supply NA here if the graph 
#' has a weight edge attribute, but you want to ignore it. Larger edge weights
#' correspond to stronger connections. This argument is not settable for 
#' "infomap" method.
#' @param steps The number of steps to take, this is actually the number of 
#' tries to make a step. It is not a particularly useful parameter. This 
#' argument is settable only for "leadingEigen"and"walktrap" method.
#' @param spins Integer constant, the number of spins to use. This is the upper 
#' limit for the number of communities. It is not a problem to supply a 
#' (reasonably) big number here, in which case some spin states will be 
#' unpopulated. This argument is settable only for "spinglass" method.
#' @param e.weights If not NULL, then a numeric vector of edge weights. 
#' The length must match the number of edges in the graph. By default the 
#' ‘weight’ edge attribute is used as weights. If it is not present, then all
#' edges are considered to have the same weight. Larger edge weights correspond 
#' to stronger connections.  This argument is settable only for "infomap"
#'  method.
#' @param v.weights If not NULL, then a numeric vector of vertex weights. The
#' length must match the number of vertices in the graph. By default the 
#' ‘weight’ vertex attribute is used as weights. If it is not present, then all
#' vertices are considered to have the same weight. A larger vertex weight means
#' a larger probability that the random surfer jumps to that vertex. This 
#' argument is settable only for "infomap" method.
#' @param nb.trials The number of attempts to partition the network (can be any
#' integer value equal or larger than 1). This argument is settable only for
#' "infomap" method.
#' @param directed Logical constant, whether to calculate directed edge 
#' betweenness for directed graphs. This argument is settable only for 
#' "edgeBetweenness" method.
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

membershipCommunities<- function(graph,
                                 method=c("walktrap", "edgeBetweenness", 
                                        "fastGreedy", "louvain", "spinglass", 
                                        "leadingEigen", "labelProp", "infomap",
                                        "optimal", "other"),
                                 FUN=NULL, directed=FALSE, weights=NULL, 
                                 steps=4, spins=25, e.weights=NULL, 
                                 v.weights=NULL, nb.trials=10)
{
    method <- match.arg(method)
    members <- membership(methodCommunity(graph=graph, method=method,
                                            FUN=FUN,
                                            directed=directed,
                                            weights=weights, 
                                            steps=steps, 
                                            spins=spins, 
                                            e.weights=e.weights, 
                                            v.weights=v.weights, 
                                            nb.trials=nb.trials))
    
    return(members)
}


################ PLOT GRAPH ###############
#' plotGraph
#'
#' @description Graphical interactive representation of the network.
#' @param graph The output of prepGraph.
#'
#' @return Creates an interactive plot, a D3 JavaScript network graph.
#' @import networkD3
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' plotGraph (graph)
plotGraph <- function(graph)
{
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph)
    plot <- networkD3::simpleNetwork(graph_d3$links, opacity=0.8, zoom=TRUE,
                                linkColour = "B1AEAE", nodeColour = "#2E66AC",
                                fontSize=12)
    return(plot)
}   


######################## PLOT COMMUNITIES ##############
#' plotComm
#' 
#' @description Graphical interactive representation of the network and its 
#' communities.
#' 
#' @param graph The output of prepGraph.
#' @param members A membership vector of the community structure, the output of
#' membershipCommunities. 
#'
#' @return Creates an interactive plot with colorful communities, a D3 
#' JavaScript network graph.
#' @import networkD3 
#' @importFrom methods is
#' @export
#'
#' @examples
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' members <- membershipCommunities (graph=graph, method="louvain")
#' plotComm(graph, members)
plotComm <- function(graph, members) 
{
    stopifnot(methods::is(graph, "igraph"))
    stopifnot(methods::is(members, "membership"))
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph, group=members)
    # Create network plot
    plot <- networkD3::forceNetwork(Links=graph_d3$links, Nodes=graph_d3$nodes,
                                    Source='source', Target='target', 
                                    NodeID='name', Group='group', opacity=0.8,
                                    fontSize=12, legend=TRUE)
    return(plot)
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
#' @param FUN see \code{\link{methodCommunity}}.
#' @param measure The measure for the comparison of the communities "vi", "nmi",
#' "split.join", "adjusted.rand"
#' @param weights this argument is not settable for "infomap" method
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" 
#' method
#' @param spins This argument is settable only for "infomap" method
#' @param e.weights This argument is settable only for "infomap" method
#' @param v.weights This argument is settable only for "infomap" method
#' @param nb.trials This argument is settable only for "infomap" method
#' @param directed This argument is settable only for "edgeBetweenness" method
#' @keywords internal

rewireCompl <- function(data, number, community, 
                        method=c("walktrap", "edgeBetweenness", 
                                 "fastGreedy", "louvain", "spinglass", 
                                 "leadingEigen", "labelProp", "infomap",
                                 "optimal", "other"),
                        FUN=NULL,
                        measure= c("vi", "nmi","split.join", "adjusted.rand"),
                        directed=FALSE, weights=NULL, steps=4, spins=25, 
                        e.weights=NULL, v.weights=NULL, nb.trials=10)
{
    method <- match.arg(method)
    measure <- match.arg (measure)
    graphRewire <- igraph::rewire(data,
                                  with=keeping_degseq(loops=FALSE,niter=number))
    comR <- membershipCommunities(graph=graphRewire, method=method, FUN=FUN,
                            directed=directed, weights=weights, steps=steps, 
                            spins=spins, e.weights=e.weights, 
                            v.weights=v.weights, nb.trials=nb.trials)
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
#' robinRobust
#'
#' @description This functions implements a procedure to examine the stability 
#' of the partition recovered by some algorithm against random perturbations 
#' of the original graph structure.
#' @param graph The output of prepGraph.
#' @param graphRandom The output of random function.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal".
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), and has to return a community object.
#' @param measure The stability measure, one of "vi", "nmi", "split.join", 
#' "adjusted.rand".
#' @param type The type of robin costruction, dependent or independent data
#' @param weights this argument is not settable for "infomap" method.
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" 
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method.
#' @param verbose flag for verbose output (default as FALSE).
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
#' measure="vi",type="independent")
robinRobust <- function(graph, graphRandom, 
                method=c("walktrap", "edgeBetweenness", 
                         "fastGreedy", "louvain", "spinglass", 
                         "leadingEigen", "labelProp", "infomap",
                         "optimal", "other"),
                FUN=NULL, measure= c("vi", "nmi","split.join", "adjusted.rand"),
                type=c("independent","dependent"), directed=FALSE, weights=NULL, 
                steps=4, spins=25, e.weights=NULL, v.weights=NULL, 
                nb.trials=10, verbose=FALSE) 
{   
    measure <- match.arg(measure)
    type<- match.arg(type)
    method <- match.arg(method)
    nrep <- 10
    comReal <- membershipCommunities(graph=graph, method=method, 
                                    FUN=FUN,
                                    directed=directed,
                                    weights=weights, 
                                    steps=steps, 
                                    spins=spins, 
                                    e.weights=e.weights, 
                                    v.weights=v.weights, 
                                    nb.trials=nb.trials) # real network

    comRandom <- membershipCommunities(graph=graphRandom, method=method,
                                        FUN=FUN,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials) # random network
    #stopifnot(length(table(comRandom))>1)
    #if(length(table(comRandom))==1) {stop("Not random graph")}
    de <- igraph::gsize(graph)
    Measure <- NULL
    vector <- NULL
    vectRandom <- NULL
    graphRewireRandom <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
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
                                    measure=measure,
                                    directed=directed,
                                    weights=weights,
                                    steps=steps, 
                                    spins=spins, 
                                    e.weights=e.weights, 
                                    v.weights=v.weights, 
                                    nb.trials=nb.trials)
                if((measure=="vi")|(measure=="split.join"))
                {
                    vector[k] <- Real$Measure
                } 
                else {
                    vector[k] <- 1-(Real$Measure)
                }
                measureReal[count2, count] <- vector[k]
                graphRewire <- Real$graphRewire
                
                #RANDOM
                Random <- rewireCompl(data=graphRandom,
                                        number=z, 
                                        community=comRandom, 
                                        method=method,
                                        measure=measure,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials)
                if((measure=="vi")|(measure=="split.join"))
                {
                    vectRandom[k] <- Random$Measure
                } else {
                    vectRandom[k] <- 1-(Random$Measure)
                }
                measureRandom[count2, count] <- vectRandom[k]
                graphRewireRandom <- Random$graphRewire
                for(k in c(2:nrep))
                {
                    count2 <- count2+1
                    Real <- rewireCompl(data=graphRewire, 
                                        number=round(0.01*z), 
                                        community=comReal, 
                                        method=method,
                                        measure=measure,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials)
                    if((measure=="vi")|(measure=="split.join"))
                    {
                        vector[k] <- Real$Measure
                    } else {
                        vector[k] <- 1-(Real$Measure)
                    }
                    measureReal[count2, count] <- vector[k]
                    Random <- rewireCompl(data=graphRewireRandom, 
                                          number=round(0.01*z),
                                          community=comRandom,
                                          method=method,
                                          measure=measure,
                                          directed=directed,
                                          weights=weights,
                                          steps=steps, 
                                          spins=spins, 
                                          e.weights=e.weights, 
                                          v.weights=v.weights, 
                                          nb.trials=nb.trials)
                    if((measure=="vi")|(measure=="split.join"))
                    {
                        vectRandom[k] <- Random$Measure
                    } else {
                        vectRandom[k] <- 1-(Random$Measure)
                    }
                    measureRandom[count2, count] <- vectRandom[k]
                }
                MeanRandom[s, count] <- mean(vectRandom)
                Mean[s, count] <- mean(vector)
            }
            if(verbose) cat("Perturbed ", z, " nodes\n")
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
                                        method=method,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials,
                                        FUN=FUN)
                Measure <- igraph::compare(comReal, comr, method=measure)
                if((measure=="vi")|(measure=="split.join"))
                {
                    vector[k] <- Measure
                } else {
                    vector[k] <- 1-Measure
                }
                measureReal1[count2] <- vector[k]
                diff <- igraph::difference(graph, graphRewire)
                
                ###RANDOM
                graphRewireRandom <- rewireOnl(data=graphRandom, number=z)
                graphRewireRandom <- igraph::union(graphRewireRandom, diffR)
                comr <- membershipCommunities(graph=graphRewireRandom,
                                        method=method,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials,
                                        FUN=FUN)
                Measure <- igraph::compare(comRandom, comr,
                                           method=measure)
                if((measure=="vi")|(measure=="split.join"))
                {
                    vectRandom[k] <- Measure
                } else {
                    vectRandom[k] <- 1-Measure
                }
                
                measureRandom1[count2] <- vectRandom[k]
                diffR <- igraph::difference(graphRandom, graphRewireRandom)
                
                for(k in c(2:nrep)) 
                {
                    count2 <- count2+1
                    ##REAL
                    Real <- rewireCompl(data=graphRewire, number=round(0.01*z),
                                        method=method,
                                        measure=measure,
                                        community=comReal,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials)
                    if((measure=="vi")|(measure=="split.join"))
                    {
                        vector[k] <- Real$Measure
                    } else {
                        vector[k] <- 1-(Real$Measure)
                    }
                    measureReal1[count2] <- vector[k]
                    
                    ## RANDOM
                    Random <- rewireCompl(data=graphRewireRandom,
                                            method=method,
                                            measure=measure,
                                            number=round(0.01*z),
                                            community=comRandom,
                                            directed=directed,
                                            weights=weights,
                                            steps=steps, 
                                            spins=spins, 
                                            e.weights=e.weights, 
                                            v.weights=v.weights, 
                                            nb.trials=nb.trials)
                    if((measure=="vi")|(measure=="split.join"))
                    {
                        vectRandom[k] <- Random$Measure
                    } else {
                        vectRandom[k] <- 1-(Random$Measure)
                    }
                    measureRandom1[count2] <- vectRandom[k]
                }
                Mean1[s] <- mean(measureReal1)
                MeanRandom1[s] <- mean(measureRandom1)  
            }
            graph <- igraph::intersection(graph, graphRewire)
            graphRandom <- igraph::intersection(graphRandom, graphRewireRandom)
            measureRandom <- cbind(measureRandom,measureRandom1)
            measureReal <- cbind(measureReal,measureReal1)
            Mean <- cbind(Mean,Mean1)
            MeanRandom <- cbind(MeanRandom,MeanRandom1)
            z1 <- igraph::gsize(graph)
            #print(z1)
            if(verbose) cat("Perturbed ", z, " nodes\n")
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



############PLOT##############
#' plotRobin
#'
#' @description This function plots two curves: the measure of the null model and the measure
#' of the real graph or the measure of the two community detection algorithms.
#' @param graph The output of prepGraph
#' @param model1 The Mean output of the robinRobust function or the Mean1 
#' output of robinCompare.
#' @param model2 The MeanRandom output of the robinRobust function or the Mean2 
#' output of robinCompare.
#' @param measure The stability measure: one of "vi", "nmi", "split.join", 
#' "adjusted.rand".
#' @param legend The legend for the graph. The default is c("model1", 
#' "model2"). If using robinRobust is recommended c("real data", "null model"), 
#' if using robinCompare, enter the names of the community detection 
#' algorithms.
#' @param title The title for the graph. The default is "Robin plot".
#'
#' @return A ggplot object.
#' @import ggplot2 igraph
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' Proc <- robinRobust(graph=graph, graphRandom=graphRandom, method="louvain",
#' type="independent")
#' plotRobin(graph=graph, model1=Proc$Mean, model2=Proc$MeanRandom,
#' measure="vi", legend=c("real data", "null model"))
#' 
plotRobin <- function(graph, model1, model2,
                      measure= c("vi", "nmi","split.join", "adjusted.rand"),
                      legend=c("model1", "model2"),
                      title="Robin plot")
{   
    measure <- match.arg (measure)
    if(measure=="vi")
    {
      N <- igraph::vcount(graph)
      mvimodel1 <- cbind(as.vector((apply(model1, 2, mean))/log2(N)),legend[1])
      mvimodel2 <- cbind(as.vector((apply(model2, 2, mean))/log2(N)),legend[2]) 
    } else if(measure=="split.join") {
      N <- igraph::vcount(graph)
      mvimodel1 <- cbind(as.vector((apply(model1, 2, mean))/(2*N)),legend[1])
      mvimodel2 <- cbind(as.vector((apply(model2, 2, mean))/(2*N)),legend[2])     
    } else {
      mvimodel1 <- cbind(as.vector((apply(model1, 2, mean))), legend[1])
      mvimodel2 <- cbind(as.vector((apply(model2, 2, mean))), legend[2])
    }
    
    percPert <- rep((seq(0,60,5)/100), 2)
    model <- mvimodel2
    mvi <- rbind(mvimodel1, model)
    colnames(mvi) <- c("mvi", "model")
    dataFrame <- data.frame(percPert, mvi)
    plot <- ggplot2::ggplot(dataFrame, aes(x=percPert, 
                                            y=as.numeric(as.character(mvi)), 
                                            colour=model, group=factor(model)))+ 
        geom_line()+
        geom_point()+ 
        xlab("Percentage of perturbation") +
        ylab("Measure") +
        ggtitle(title)
      
    return(plot)
}

############### COMPARISON DIFFERENT METHODS ##########
#' robinCompare
#'
#' @description  This function compares the robustness of two community 
#' detection algorithms.
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
#' "adjusted.rand".
#' @param type The type of robin costruction, dependent or independent.
#' @param weights This argument is not settable for "infomap" method.
#' @param steps This argument is settable only for "leadingEigen"and"walktrap" 
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method.
#' @param verbose flag for verbose output (default as FALSE).
#' 
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method 
#' - the matrix "Mean2" with the means of the procedure for the second method.
#' 
#' @import igraph
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' robinCompare(graph=graph, method1="louvain", 
#' method2="fastGreedy", measure="vi", type="independent")
robinCompare <- function(graph, 
                      method1=c("walktrap", "edgeBetweenness", "fastGreedy",
                                "leadingEigen","louvain","spinglass",
                                "labelProp","infomap","optimal", "other"),
                      method2=c("walktrap", "edgeBetweenness", "fastGreedy",
                                "leadingEigen","louvain","spinglass",
                                "labelProp","infomap","optimal", "other"),
                      FUN1=NULL, FUN2=NULL,
                      measure= c("vi", "nmi","split.join", "adjusted.rand"),
                      type=c("independent", "dependent"),
                      directed=FALSE, weights=NULL, steps=4, 
                      spins=25, e.weights=NULL, v.weights=NULL, 
                      nb.trials=10, verbose=FALSE)
{   
    method1 <- match.arg(method1)
    method2 <- match.arg(method2)
    type <- match.arg(type)
    measure <-match.arg(measure)
    nrep <- 10
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
    Measure <- NULL
    vector1 <- NULL
    vector2 <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detected robin method ", type, " type\n")
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
                comr1 <- membershipCommunities(graph=graphRewire, 
                                               method=method1,
                                               FUN=FUN1,
                                               directed=directed,
                                               weights=weights,
                                               steps=steps, 
                                               spins=spins, 
                                               e.weights=e.weights, 
                                               v.weights=v.weights, 
                                               nb.trials=nb.trials)
                comr2 <- membershipCommunities(graph=graphRewire, 
                                               method=method2, 
                                               FUN=FUN2,
                                               directed=directed,
                                               weights=weights,
                                               steps=steps, 
                                               spins=spins, 
                                               e.weights=e.weights, 
                                               v.weights=v.weights, 
                                               nb.trials=nb.trials)
                if((measure=="vi")|(measure=="split.join"))
                {
                    vector1[k] <- igraph::compare(comr1, comReal1, 
                                                  method=measure)
                    vector2[k] <- igraph::compare(comr2, comReal2, 
                                                  method=measure)
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
                    comr1 <- membershipCommunities(graph=graphRewire,
                                                   method=method1,
                                                   FUN=FUN1,
                                                   directed=directed,
                                                   weights=weights,
                                                   steps=steps, 
                                                   spins=spins, 
                                                   e.weights=e.weights, 
                                                   v.weights=v.weights, 
                                                   nb.trials=nb.trials)
                    comr2 <- membershipCommunities(graph=graphRewire, 
                                                   FUN=FUN2,
                                                   method=method2,
                                                   directed=directed,
                                                   weights=weights,
                                                   steps=steps, 
                                                   spins=spins, 
                                                   e.weights=e.weights, 
                                                   v.weights=v.weights, 
                                                   nb.trials=nb.trials)
                    if((measure=="vi")|(measure=="split.join"))
                    {
                        vector1[k] <- igraph::compare(comr1, comReal1, 
                                                      method=measure)
                        vector2[k] <- igraph::compare(comr2, comReal2, 
                                                      method=measure)
                    } else {
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
            if(verbose) cat("Perturbed ", z, " nodes\n")
        }
        
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
                comr1 <- membershipCommunities(graph=graphRewire, 
                                               method=method1,
                                               FUN=FUN1,
                                               directed=directed,
                                               weights=weights,
                                               steps=steps, 
                                               spins=spins, 
                                               e.weights=e.weights, 
                                               v.weights=v.weights, 
                                               nb.trials=nb.trials)
                comr2 <- membershipCommunities(graph=graphRewire, 
                                               method=method2,
                                               FUN=FUN2,
                                               directed=directed,
                                               weights=weights,
                                               steps=steps, 
                                               spins=spins, 
                                               e.weights=e.weights, 
                                               v.weights=v.weights, 
                                               nb.trials=nb.trials)
                if((measure=="vi")|(measure=="split.join"))
                {
                    vector1[k] <- igraph::compare(comr1, comReal1, 
                                                  method= measure)
                    vector2[k] <- igraph::compare(comr2, comReal2, 
                                                  method= measure)
                } else {
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
                    comr1 <- membershipCommunities(graph=graphRewire,
                                                   method=method1,
                                                   FUN=FUN1,
                                                   directed=directed,
                                                   weights=weights,
                                                   steps=steps, 
                                                   spins=spins, 
                                                   e.weights=e.weights, 
                                                   v.weights=v.weights, 
                                                   nb.trials=nb.trials)
                    comr2 <- membershipCommunities(graph=graphRewire,
                                                   method=method2,
                                                   FUN=FUN2,
                                                   directed=directed,
                                                   weights=weights,
                                                   steps=steps, 
                                                   spins=spins, 
                                                   e.weights=e.weights, 
                                                   v.weights=v.weights, 
                                                   nb.trials=nb.trials)
                    if((measure=="vi")|(measure=="split.join"))
                    {
                        vector1[k] <- igraph::compare(comr1, comReal1, 
                                                      method=measure)
                        vector2[k] <- igraph::compare(comr2, comReal2, 
                                                      method=measure)
                    } else {
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
            if(verbose) cat("Perturbed ", z, " nodes\n")
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
#' @param measure The measure for the comparison of the communities "vi", "nmi",
#' "split.join", "adjusted.rand".
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
                                measure= c("vi", "nmi","split.join", 
                                           "adjusted.rand"),
                                muParam=0, orderParam=4, nKnots=7, 
                                BParam=10000, isPaired=TRUE) 
{
    measure <- match.arg (measure)
    if(measure=="vi")
    {
    nOrder <- log2(igraph::gorder(graph))
    modeled1 <- as.matrix(model1)/nOrder
    modeled2 <- as.matrix(model2)/nOrder
    }else if(measure=="split.join")
    {N <- igraph::vcount(graph)
     modeled1 <- as.matrix(model1)/(2*N)
     modeled2 <- as.matrix(model2)/(2*N)
    }else{
    modeled1 <- as.matrix(model1)
    modeled2 <- as.matrix(model2)   
    }
   
    
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
#' @param model1 The Mean output of the robinRobust function (or the Mean1 
#' output of the robinCompare function).
#' @param model2 The MeanRandom output of the robinRobust function (or the 
#' Mean2 output of the robinCompare function).
#' @param verbose flag for verbose output (default as FALSE).
#' 
#' @return A numeric value, the Bayes factor
#' @import gprege 
#' @importFrom stats sd var
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' Proc <- robinRobust(graph=graph, graphRandom=graphRandom, 
#' method="louvain", measure="vi",type="independent")
#' robinGPTest(model1=Proc$Mean,model2=Proc$MeanRandom)
robinGPTest <- function(model1, model2, verbose=FALSE)
{ 
   ratios <- log2((model1+0.001)/(model2+0.001))
   #rapporto tra la media delle misure tra il modello reale e quello perturbato 
   #e la media delle distanze tra il random e la sua perturbazione
   res <- as.vector(ratios)

   nRewire <- seq(0,60,5)
   nrep <- 10
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
    gpregeOutput <- gprege::gprege(data=datadum, inputs=dvet,
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
#' @param graph The output of prepGraph.
#' @param model1 The Mean output of the robinRobust function (or the Mean1 
#' output of the robinCompare function).
#' @param model2 The MeanRandom output of the robinRobust function (or the 
#' Mean2 output of the robinCompare function).
#' @param measure The stability measure "vi", "nmi", "split.join", 
#' "adjusted.rand".
#' @param legend The legend for the graph. The default is c("model1", 
#' "model2"). If using robinRobust is recommended c("real data", "null model"), 
#' if using robinCompare, enter the names of the community detection 
#' algorithms.
#' @param verbose flag for verbose output (default as FALSE).
#' 
#' @return Two plots: the fitted curves and the adjusted p-values. A vector of the adjusted p-values. 
#' @import igraph ggplot2 fdatest graphics qpdf
#' @export
#'
#' @examples 
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' Proc <- robinRobust(graph=graph, graphRandom=graphRandom, method="louvain",
#' measure="vi",type="independent")
#' robinFDATest(graph=graph, model1=Proc$Mean, model2=Proc$MeanRandom, 
#' measure="vi",legend=c("real data", "null model"))
robinFDATest <- function(graph,model1,model2, measure= c("vi", "nmi",
                        "split.join", "adjusted.rand"),
                        legend=c("model1", "model2"), verbose=FALSE)
{
    if(verbose) cat("Computing Interval testing procedure.\n")
    object <- createITPSplineResult(graph, model1, model2, measure)
    
    
  #Functional Data plot
    
    J <- dim(object$data.eval)[2]
    xmin <- 0
    xmax <- 0.6
    Abscissa <- seq(xmin,xmax,len=J)
    
    model1 <- cbind(as.numeric(as.vector(t(object$data.eval[1:10,]))), 
                    legend[1], rep(1:10,each=1000), as.numeric(rep(Abscissa, 
                                                                   times = 10)))
    model2 <- cbind(as.numeric(as.vector(t(object$data.eval[11:20,]))), 
                    legend[2], rep(11:20,each=1000), as.numeric(rep(Abscissa, 
                                                                    times = 10))
                    )
    measure <- rbind(model1, model2)
    colnames(measure) <- c("measure","model","s","percPert")
    dataFrame <- data.frame(measure)
    plot1<-ggplot2::ggplot(dataFrame, ggplot2::aes(x=as.numeric(percPert),
                 y=as.numeric(measure), color= model, group=s)) +
              ggplot2::geom_line() +
              ggplot2::xlab("Percentage of perturbation") +
              ggplot2::ylab("Measure")+
              ggplot2::ggtitle("Functional Data Analysis")+
              ggplot2::scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 
                                                     0.5, 0.6))  
    
     
  
    #P value plot
     p <- length(object$pval)
     xmin <- 0
     xmax <- 0.6
     abscissa.pval <- rep(seq(xmin,xmax,len=p),time=2)
     pvalue <- c(object$pval,object$corrected.pval)
     type <- c(rep("pvalue",p),rep("pvalue.adj",p))
     PdataFrame<-data.frame(cbind(abscissa.pval,pvalue,type))
     plot2<-ggplot2::ggplot(PdataFrame, ggplot2::aes(x=as.numeric(abscissa.pval),
                                                    y=as.numeric(pvalue), color= type)) +
       ggplot2::geom_point() +
       ggplot2::xlab("Percentage of perturbation") +
       ggplot2::ylab("p_value")+
       ggplot2::ggtitle("P-values")+
       ggplot2::scale_x_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 
                                              0.5, 0.6))+
       ggplot2::geom_hline(yintercept = 0.05,color = "red")
     
    
    plot <- gridExtra::grid.arrange(plot1,plot2, ncol=2)
    
    print(plot)
   
    adj.pvalue<-object$corrected.pval
    pvalue<-object$pval
    output<-list(adj.pvalue=adj.pvalue,
                 pvalues=pvalue)
    return(output)
}  

########### AREA UNDER THE CURVE    ##############
#' robinAUC
#'
#' @description This function calculates the area under two curves with a spline approach. 
#' @param graph The output of prepGraph.
#' @param model1 The Mean output of the robinRobust function (or the Mean1 
#' output of the robinCompare function).
#' @param model2 The MeanRandom output of the robinRobust function (or the 
#' Mean2 output of the robinCompare function).
#' @param measure The stability measure "vi", "nmi", "split.join", 
#' "adjusted.rand".
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
#' Proc <- robinRobust(graph=graph, graphRandom=graphRandom, method="louvain",
#' measure="vi",type="independent")
#' robinAUC(graph=graph, model1=Proc$Mean, model2=Proc$MeanRandom)
robinAUC <- function(graph, model1, model2, 
                        measure= c("vi", "nmi","split.join", "adjusted.rand"),
                        verbose=FALSE)
{
    if(verbose) cat("Computing area under the curve (AUC).\n")
    measure <- match.arg (measure)
    if(measure=="vi")
    {
        N <- igraph::vcount(graph)
        mvimeanmodel1 <- cbind(as.vector((apply(model1, 2, mean))/log2(N)))
        mvimeanmodel2 <- cbind(as.vector((apply(model2, 2, mean))/log2(N))) 
    }else if(measure=="split.join")
    {
        N <- igraph::vcount(graph)
        mvimeanmodel1 <- cbind(as.vector((apply(model1, 2, mean))/(2*N)))
        mvimeanmodel2 <- cbind(as.vector((apply(model2, 2, mean))/(2*N)))     
    }else
    {
        mvimeanmodel1 <- cbind(as.vector((apply(model1, 2, mean))))
        mvimeanmodel2 <- cbind(as.vector((apply(model2, 2, mean))))
    }
    area1 <- DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel1, 
                          method ="spline")
    area2 <- DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel2, 
                          method ="spline")
    output <- list(area1=area1,area2=area2)
return(output)
}

