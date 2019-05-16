######PREPARATION GRAPH########## 

#' prepGraph
#' 
#' @description The prepGraph function is able to read graphs from a file and 
#' to prepare them for the analysis.
#'
#' @param file The file to read from.
#' @param file.format Character constant giving the file format. Right now
#' as_edgelist, pajek, graphml, gml, ncol, lgl, dimacs, graphdb and igraph are
#' supported
#' @param numbers A logical value indicating if the names of the nodes are 
#' values.This argument is settable for the edgelist format. 
#' The default is FALSE.
#' @param header A logical value indicating whether the file contains 
#' the names of the variables as its first line.This argument is settable 
#' for the edgelist format.The default is FALSE.
#' @return A graph which do not contain loop and multiple edges.
#' @import igraph
#' @export
#'
#' @examples 
#' graph <- prepGraph(file=my_file, file.format="gml")
prepGraph <- function(file,
                        file.format=c("edgelist", "pajek", "ncol", "lgl",
                                      "graphml", "dimacs", "graphdb", "gml",
                                      "dl","igraph"),
                        numbers= FALSE,
                        directed=FALSE,
                        header=FALSE)
{ 
    file.format <- match.arg(file.format)
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
        edge <- read.table(file,colClasses = "character", quote="\"",
                           header=header)
        edge <- as.matrix(edge)
        graph <- igraph::graph_from_edgelist(edge, directed=directed)
        graph <- igraph::simplify(graph)
    }else
    {
        net <- igraph::read_graph(file=file, format=file.format,
                                  directed=directed
                                  )
        ind <- igraph::V(net)[degree(net) == 0] #isolate node
        graph <- igraph::delete.vertices(net, ind)
        graph <- igraph::simplify(graph)
    }
    ##Other method:
    # method <- match.arg(method)
    # if(method == "igraph")
    # {
    #     graph <- igraph::read_graph(file, format=file.format, 
    #                                 directed=is.directed)
    # } else if(method == "robin") {
    #     edge <- read.table(file, quote="\"")
    #     edge <- as.matrix(edge)
    #     
    #     #graph_from_edgelist(edge,direct=FALSE)
    #     
    #     vet1 <- as.vector(t(edge))
    #     un <- unique(sort(vet1))
    #     if(max(un) != length(un)) 
    #     { ##se non vi sono tutti i nodi il massimo del vettore 
    #         #non è uguale alla lunghezza
    #         id <- seq(1, length(un))  ##il nome dei nodi è vet 
    #         #che è uguale all'id, altrimenti è vet1
    #         vet <- vet1
    #         for(i in c(1:length(un))) 
    #         {
    #             ind <- which(vet1 == un[i])
    #             vet[ind] <- id[i]
    #         }
    #         edge <- matrix(vet, ncol=2, byrow=TRUE)
    #         graph <- igraph::graph(vet, directed=is.directed)
    #     } else {
    #         graph <- igraph::graph(vet1, directed=is.directed)
    #     }
    # }
    # graph <- igraph::simplify(graph) 
    return(graph)
}


######GRAPH RANDOM#########
#' random
#'
#' @description Randomly rewire the edges while preserving the original graph's 
#' degree distribution.
#' @param graph The output of prepGraph.
#'
#' @return A randomly rewired graph.
#' @import igraph
#' @export
#'
#' @examples 
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
random <- function(graph)
{
    z <- igraph::gsize(graph) ## number of edges
    graphRandom <- igraph::rewire(graph, 
                            with=igraph::keeping_degseq(loops=FALSE, niter=z))
    #rewiring for z all the edges
    return(graphRandom)
}


#####COMMUNITY METHOD####    
#' methodCommunity
#' 
#' @description This function detects the community structure of a graph.
#' The community structure are found by all functions implemented in igraph.
#' @param graph The output of prepGraph.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal", "other".
#' @param FUN see \code{\link{methodCommunity}}.
#' @param FUN in case the @method parameter is "other" there is the possibility 
#' to use a personal function passing its name through this parameter.
#' The personal parameter has to take as input the @graph and the @weights 
#' (that can be NULL), moreover it has to return a communities object.
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
#' @return a Communities object.
#' @import igraph
#' @export
#'
#' @examples 
#' graph <- prepGraph(file=my_file, file.format="gml")
#' methodCommunity (graph=graph, method="louvain") 
methodCommunity <- function(graph, 
                            method=c("walktrap", "edgeBetweenness", 
                                    "fastGreedy", "louvain", "spinglass", 
                                    "leadingEigen", "labelProp", "infomap",
                                    "optimal", "other"),
                            FUN=NULL,
                            directed=FALSE,
                            weights=NULL, 
                            steps=4, 
                            spins=25, 
                            e.weights=NULL, 
                            v.weights=NULL, 
                            nb.trials=10)
{
    method <- match.arg(method)
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
            optimal=igraph::cluster_optimal(graph, 
                                            weights = weights),
            louvain=igraph::cluster_louvain(graph=graph, 
                                    weights=weights), 
           
            walktrap=igraph::cluster_walktrap(graph=graph, 
                                    weights=weights, 
                                    steps=steps), 
           
            spinglass=igraph::cluster_spinglass(graph=graph, 
                                        weights=weights, 
                                        spins=spins), 
           
            leadingEigen=igraph::cluster_leading_eigen(graph=graph, 
                                            steps=steps, 
                                            weights=weights,
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

#####MEMBERSHIP COMMUNITIES####    
#' membershipCommunities
#' 
#' @description This function gives the membership vector of the community 
#' structure. The community structure are found by all functions implemented 
#' in igraph.
#' @param graph The output of prepGraph.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal".
#' @param FUN see \code{\link{methodCommunity}}.
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
#' @return Membership vector of the community structure.
#' @import igraph
#' @export
#'
#' @examples 
#' graph <- prepGraph(file=my_file, file.format="gml")
#' membershipCommunities (graph=graph, method="louvain")

membershipCommunities<- function(graph,
                                 method=c("walktrap", "edgeBetweenness", 
                                        "fastGreedy", "louvain", "spinglass", 
                                        "leadingEigen", "labelProp", "infomap",
                                        "optimal", "other"),
                                 FUN=NULL,
                                 directed=FALSE,
                                 weights=NULL, 
                                 steps=4, 
                                 spins=25, 
                                 e.weights=NULL, 
                                 v.weights=NULL, 
                                 nb.trials=10)
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
#' @return An interactive plot.
#' @import networkD3
#' @export
#'
#' @examples 
#' graph <- prepGraph(file=my_file, file.format="gml")
#' plotGraph (graph)
plotGraph <- function(graph)
{
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph)
    plot <- networkD3::simpleNetwork(graph_d3$links, opacity=0.8, zoom=TRUE,
                                     fontSize=12)
    return(plot)
}   


######################## PLOT COMMUNITIES ##############
#' plotCommu
#' 
#' @description Graphical interactive representation of the network and its 
#' communities
#' 
#' @param graph The output of prepGraph.
#' @param members A membership vector of the community structure, the output of
#' membershipCommunities. 
#'
#' @return An interactive plot.
#' @import networkD3
#' @export
#'
#' @examples
#' graph <- prepGraph(file=my_file, file.format="gml")
#' members <- membershipCommunities (graph=graph, method="louvain")
#' plotCommu(graph, members)
plotCommu <- function(graph, members) 
{
    stopifnot(is(graph, "igraph"))
    stopifnot(is(members, "membership"))
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph, group=members)
    # Create force directed network plot
    plot <- networkD3 ::forceNetwork(Links=graph_d3$links, 
                                     Nodes=graph_d3$nodes,
                                     Source='source', 
                                     Target='target', 
                                     NodeID='name', 
                                     Group='group',
                                     opacity=0.8,
                                     fontSize=12,
                                     legend=TRUE)
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
#' @keyword internal
#' @return A list object
rewireCompl <- function(data, number, community, 
                        method=c("walktrap", "edgeBetweenness", 
                                 "fastGreedy", "louvain", "spinglass", 
                                 "leadingEigen", "labelProp", "infomap",
                                 "optimal", "other"),
                        FUN=NULL,
                        measure= c("vi", "nmi","split.join", "adjusted.rand"),
                        directed=FALSE,
                        weights=NULL, 
                        steps=4, 
                        spins=25, 
                        e.weights=NULL, 
                        v.weights=NULL, 
                        nb.trials=10)
{
    method <- match.arg(method)
    measure <- match.arg (measure)
    graphRewire <- igraph::rewire(data,
                                  with=keeping_degseq(loops=FALSE,niter=number))
    comR <- membershipCommunities(graph=graphRewire, method=method, FUN=FUN,
                            directed=directed,
                            weights=weights,
                            steps=steps, 
                            spins=spins, 
                            e.weights=e.weights, 
                            v.weights=v.weights, 
                            nb.trials=nb.trials)
    Measure <- igraph::compare(community, comR, method=measure)
    output <- list(Measure=Measure, graphRewire=graphRewire)
    return(output)
}

######### REWIRE ONLY ###########
#' rewireOnl
#' @description makes the rewire function of igraph
#' @param data The output of prepGraph
#' @param number Number of rewiring trials to perform.
#' @keyword internal
#' @return An igraph object
rewireOnl <- function(data, number)
{
    graphRewire <- igraph::rewire(data,
                                    with=keeping_degseq(loops=FALSE,
                                                        niter=number))
    return(graphRewire)
}


########  ROBIN PROCEDURE #######
#' robinProc
#'
#' @description A procedure to examine the stability of the partition recovered 
#' against random perturbations of the original graph structure.
#' @param graph The output of prepGraph.
#' @param graphRandom The output of random function.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal".
#' @param FUN see \code{\link{methodCommunity}}.
#' @param measure The measure for the comparison of the communities "vi", "nmi",
#' "split.join", "adjusted.rand"
#' @param type The type of robin costruction dependent or independent data
#' @param weights this argument is not settable for "infomap" method.
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" 
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method.
#'
#' @return A list object.
#' @import igraph
#' @export
#'
#' @examples 
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' robinProc(graph=graph, graphRandom=graphRandom, method="louvain",
#' measure="vi",type="independent")
robinProc <- function(graph, graphRandom, 
                method=c("walktrap", "edgeBetweenness", 
                         "fastGreedy", "louvain", "spinglass", 
                         "leadingEigen", "labelProp", "infomap",
                         "optimal", "other"),
                FUN=NULL,
                measure= c("vi", "nmi","split.join", "adjusted.rand"),
                type=c("dependent", "independent"),
                directed=FALSE,
                weights=NULL, 
                steps=4, 
                spins=25, 
                e.weights=NULL, 
                v.weights=NULL, 
                nb.trials=10) 
{
    measure <- match.arg (measure)
    method <- match.arg(method)
    type <- match.arg(type)
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
    de <- igraph::gsize(graph)
    Measure <- NULL
    vector <- NULL
    vectRandom <- NULL
    graphRewireRandom <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
    #INDEPENDENT    
    if(type == "independent") 
    {
        #OUTPUT MATRIX
        measReal <- matrix(0, nrep^2, length(nRewire))
        Random <- matrix(0, nrep^2, length(nRewire))
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
                                    directed=directed,
                                    weights=weights,
                                    steps=steps, 
                                    spins=spins, 
                                    e.weights=e.weights, 
                                    v.weights=v.weights, 
                                    nb.trials=nb.trials)
                if((measure=="vi")|(measure=="split.join"))
                {
                    vector[k] <- Real$Meaure
                } else {
                    vector[k] <- 1-(Real$Meaure)
                }
                measReal[count2, count] <- vector[k]
                graphRewire <- Real$graphRewire
                
                #RANDOM
                Random <- rewireCompl(data=graphRandom,
                                        number=z, 
                                        community=comRandom, 
                                        method=method, 
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
                Random[count2, count] <- vectRandom[k]
                graphRewireRandom <- Random$graphRewire
                for(k in c(2:nrep))
                {
                    count2 <- count2+1
                    Real <- rewireCompl(data=graphRewire, 
                                        number=round(0.01*z), 
                                        community=comReal, 
                                        method=method,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials)
                    if((measure=="vi")|(measure=="split.join"))
                    {
                        vector[k] <- Real$Meaure
                    } else {
                        vector[k] <- 1-(Real$Meaure)
                    }
                    measReal[count2, count] <- vector[k]
                    Random <- rewireCompl(data=graphRewireRandom, 
                                          number=round(0.01*z),
                                          community=comRandom,
                                          method=method,
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
                    Random[count2, count] <- vectRandom[k]
                }
                MeanRandom[s, count] <- mean(vectRandom)
                Mean[s, count] <- mean(vector)
            }
            print(z) 
        }
  #DEPENDENT 
    }else{
        z <- round((5*de)/100, 0) #the 5% of the edges
        measReal <- rep(0, nrep^2)
        measReal1 <- NULL
        Random <- rep(0, nrep^2)
        Random1 <- NULL
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
                if((measure=="vi")|(measure=="split.join"))
                {
                    vector[k] <- igraph::compare(comReal, comr, method=measure)
                } else {
                    vector[k] <- 1-(igraph::compare(comReal, comr, 
                                                    method=measure))
                }
                measReal1[count2] <- vector[k]
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
                if((measure=="vi")|(measure=="split.join"))
                {
                    vectRandom[k] <- igraph::compare(comRandom, comr, 
                                                     method=measure)
                } else {
                    vectRandom[k] <- 1-(igraph::compare(comRandom, comr, 
                                                        method=measure))
                }
                
                Random1[count2] <- vectRandom[k]
                diffR <- igraph::difference(graphRandom, graphRewireRandom)
                
                for(k in c(2:nrep)) 
                {
                    count2 <- count2+1
                    ##REAL
                    Real <- rewireCompl(data=graphRewire, number=round(0.01*z),
                                        method=method,
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
                        vector[k] <- Real$Meaure
                    } else {
                        vector[k] <- 1-(Real$Meaure)
                    }
                    measReal1[count2] <- vector[k]
                    
                    ## RANDOM
                    Random <- rewireCompl(data=graphRewireRandom,
                                            method=method,
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
                    Random1[count2] <- vectRandom[k]
                }
                Mean1[s] <- mean(measReal1)
                MeanRandom1[s] <- mean(Random1)  
            }
            graph <- igraph::intersection(graph, graphRewire)
            graphRandom <- igraph::intersection(graphRandom, graphRewireRandom)
            Random <- cbind(Random,Random1)
            measReal <- cbind(measReal,measReal1)
            Mean <- cbind(Mean,Mean1)
            MeanRandom <- cbind(MeanRandom,MeanRandom1)
            z1 <- igraph::gsize(graph)
            print(z1)
            }
    }
    colnames(Random) <- nRewire
    colnames(measReal) <- nRewire
    colnames(MeanRandom) <- nRewire
    colnames(Mean) <- nRewire
    nn <- rep(nRewire, each=nrep) 
    ratios <- log2((Mean+0.001)/(MeanRandom+0.001))
    #rapporto tra la media delle misure tra il modello reale e quello perturbato 
    #e la media delle distanze tra il random e la sua perturbazione
    bats <- as.vector(ratios)
    names(bats) <- nn
    res <- cbind(ID="ratios", t(bats))#la trasposta del rapporto
    
    output <- list( measure=measReal,
                    Random=Random,
                    Mean=Mean,
                    MeanRandom=MeanRandom,
                    ratios=res
                    )
      return(output)

}



############PLOT##############
#' plotRobin
#'
#' @description The plot of the two VI curves, the VI of the null model and 
#' of the real graph.
#' @param graph The output of prepGraph
#' @param model The viMean output of the robinProc function.
#' @param modelR The viMeanRandom output of the robinProc function.
#' @param measure The measure for the comparison of the communities "vi", "nmi",
#' "split.join", "adjusted.rand"
#' @param legend The legend for the graph. The default is c("real data", 
#' "null model").
#'
#' @return A ggplot object.
#' @import ggplot2 igraph
#' @export
#'
#' @examples 
#' graph <- prepGraph(file=my_file, file.format="gml")
#' graphRandom <- random(graph=graph)
#' Proc <- robinProc(graph=graph, graphRandom=graphRandom, method="louvain",
#' type="independent")
#' plotRobin(graph=graph, model=Proc$Mean, modelR=Proc$MeanRandom, 
#' legend=c("real data", "null model"))
plotRobin <- function(graph,
                      model,
                      modelR,
                      measure= c("vi", "nmi","split.join", "adjusted.rand"),
                      legend=c("real data", "null model"))
{   
    measure <- match.arg (measure)
    if(measure=="vi")
   {
    N <- igraph::vcount(graph)
    mvimodel1 <- cbind(as.vector((apply(model, 2, mean))/log2(N)),legend[1])
    mvimodel2 <- cbind(as.vector((apply(modelR, 2, mean))/log2(N)),legend[2]) 
    }else if(measure=="split.join")
    {
    N <- igraph::vcount(graph)
    mvimodel1 <- cbind(as.vector((apply(model, 2, mean))/(2*N)),legend[1])
    mvimodel2 <- cbind(as.vector((apply(modelR, 2, mean))/(2*N)),legend[2])     
    }else
    {
    mvimodel1 <- cbind(as.vector((apply(model, 2, mean))),legend[1])
    mvimodel2 <- cbind(as.vector((apply(modelR, 2, mean))),legend[2])
    }
    
    percPert <- rep((seq(0,60,5)/100),2)
    mvi <- rbind(mvimodel1,mvimodel2)
    colnames(mvi) <- c("mvi","model")
    dataFrame <- data.frame(percPert,mvi)
    plot <- ggplot2::ggplot(dataFrame, aes(x=percPert, 
                                            y=as.numeric(as.character(mvi)), 
                                            colour=model,
                                            group=factor(model)))+ 
        geom_line()+
        geom_point()+ 
        xlab("Percentage of perturbation") +
        ylab("Measure") +
        ggtitle("Robin plot")
        #scale_y_continuous(limits=c(0,0.6)) #only for VI
    return(plot)
}

############### COMPARISON DIFFERENT METHODS ##########
#' comparison
#'
#' @description  A procedure to compare two different methods of community 
#' detection.
#' @param graph The output of prepGraph.
#' @param graphRandom The output of random function.
#' @param method1 The first clustering method, one of "walktrap", 
#' "edgeBetweenness", "fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal".
#' @param method2 The second custering method one of "walktrap",
#' "edgeBetweenness","fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal".
#' @param FUN1 its a personal designed function when method1 is "others". 
#' see \code{\link{methodCommunity}}.
#' @param FUN2 its a personal designed function when method2 is "others". 
#' see \code{\link{methodCommunity}}.
#' @param measure The measure for the comparison of the communities "vi", "nmi",
#' "split.join", "adjusted.rand"
#' @param type The type of robin costruction dependent or independent.
#' @param weights This argument is not settable for "infomap" method.
#' @param steps This argument is settable only for "leadingEigen"and"walktrap" 
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method.
#'
#' @return A list object
#' @import igraph
#' @export
#'
#' @examples 
#' Comp <- comparison(graph=graph, graphRandom=graphRandom, 
#' method1="louvain", method2="fastGreedy",type="independent")
comparison <- function(graph,graphRandom, 
                       method1=c("walktrap", "edgeBetweenness", "fastGreedy",
                                "leadingEigen","louvain","spinglass",
                                "labelProp","infomap","optimal", "other"),
                       method2=c("walktrap", "edgeBetweenness", "fastGreedy",
                                "leadingEigen","louvain","spinglass",
                                "labelProp","infomap","optimal", "other"),
                        FUN1=NULL,
                        FUN2=NULL,
                        measure= c("vi", "nmi","split.join", "adjusted.rand"),
                        type=c("dependent", "independent"),
                        directed=FALSE,
                        weights=NULL, 
                        steps=4, 
                        spins=25, 
                        e.weights=NULL, 
                        v.weights=NULL, 
                        nb.trials=10)
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
    comRandom1 <- membershipCommunities(graph=graphRandom, method=method1, 
                                        FUN=FUN1,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials) 
    comRandom2 <- membershipCommunities(graph=graphRandom, method=method2, 
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
    vectorR1 <- NULL
    vectorR2 <- NULL
    graphRewire <- NULL
    Random<-NULL
    count <- 1
    nRewire <- seq(0,60,5)
    if(type == "independent") 
    {
        measReal1 <- matrix(0, nrep^2, length(nRewire))
        measReal2 <- matrix(0, nrep^2, length(nRewire))
        R1 <- matrix(0, nrep^2, length(nRewire))
        R2 <- matrix(0, nrep^2, length(nRewire))
        Mean1 <- matrix(0, nrep, length(nRewire))
        Mean2 <- matrix(0, nrep, length(nRewire))
        MeanRandom1 <- matrix(0, nrep, length(nRewire))
        MeanRandom2 <- matrix(0, nrep, length(nRewire))
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
                measReal1[count2, count] <- vector1[k]
                measReal2[count2, count] <- vector2[k]
                Random <- rewireOnl(data=graphRandom, number=z)
                comrR1 <- membershipCommunities(graph=Random, method=method1,
                                                FUN=FUN1,
                                                directed=directed,
                                                weights=weights,
                                                steps=steps, 
                                                spins=spins, 
                                                e.weights=e.weights, 
                                                v.weights=v.weights, 
                                                nb.trials=nb.trials)
                comrR2 <- membershipCommunities(graph=Random, method=method2,
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
                    vectorR1[k] <- igraph::compare(comrR1, comRandom1, 
                                                   method=measure)
                    vectorR2[k] <- igraph::compare(comrR2, comRandom2, 
                                                   method=measure)
                } else {
                    vectorR1[k] <- 1-(igraph::compare(comrR1, comRandom1, 
                                                   method=measure))
                    vectorR2[k] <- 1-(igraph::compare(comrR2, comRandom2, 
                                                   method=measure))
                }
                
                R1[count2, count] <- vectorR1[k]
                R2[count2, count] <- vectorR2[k]
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
                    
                    measReal1[count2, count] <- vector1[k]
                    measReal2[count2, count] <- vector2[k]
                    
                    graphRewireRandom <- rewireOnl(data=Random,
                                                    number=round(0.01*z))
                    comrR1 <- membershipCommunities(graph=Random, 
                                                    method=method1,
                                                    FUN=FUN1,
                                                    directed=directed,
                                                    weights=weights,
                                                    steps=steps, 
                                                    spins=spins, 
                                                    e.weights=e.weights, 
                                                    v.weights=v.weights, 
                                                    nb.trials=nb.trials)
                    comrR2 <- membershipCommunities(graph=Random, 
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
                        vectorR1[k] <- igraph::compare(comrR1, comRandom1, 
                                                       method=measure)
                        vectorR2[k] <- igraph::compare(comrR2, comRandom2, 
                                                       method=measure)
                    } else {
                        vectorR1[k] <- 1-(igraph::compare(comrR1, comRandom1, 
                                                       method=measure))
                        vectorR2[k] <- 1-(igraph::compare(comrR2, comRandom2, 
                                                       method=measure))
                    }
                    R1[count2, count] <- vectorR1[k]
                    R2[count2, count] <- vectorR2[k]
                }
                Mean1[s, count] <- mean(vector1)
                Mean2[s, count] <- mean(vector2)
                MeanRandom1[s, count] <- mean(vectorR1)
                MeanRandom2[s, count] <- mean(vectorR2)
            }
            print(z)
        }
        
    }else{
        z <- round((5*de)/100, 0)
        measReal1 <- rep(0, nrep^2)
        measReal11 <- NULL
        R11<-NULL
        measReal2 <- rep(0, nrep^2)
        measReal22 <- NULL
        R22 <- NULL
        Mean1 <- rep(0, nrep)
        Mean2 <-rep(0, nrep)
        Mean11 <- NULL
        Mean22 <- NULL
        MeanRandom1 <- rep(0, nrep)
        MeanRandom2 <-rep(0, nrep)
        MeanRandom11 <- NULL
        MeanRandom22 <- NULL
        diff <- NULL
        diffR <- NULL
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
             
                measReal11[count2] <- vector1[k]
                measReal22[count2] <- vector2[k]
                diff <- igraph::difference(graph, graphRewire)
                
                #Random
                Random <- rewireOnl(data=graphRandom, number=z)
                Random <- igraph::union(Random, diffR)
                comrR1 <- membershipCommunities(graph=Random, method=method1,
                                                FUN=FUN1,
                                                directed=directed,
                                                weights=weights,
                                                steps=steps, 
                                                spins=spins, 
                                                e.weights=e.weights, 
                                                v.weights=v.weights, 
                                                nb.trials=nb.trials)
                comrR2 <- membershipCommunities(graph=Random, method=method2,
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
                    vectorR1[k] <- igraph::compare(comrR1, comRandom1, 
                                                   method=measure)
                    vectorR2[k] <- igraph::compare(comrR2, comRandom2, 
                                                   method=measure)
                } else {
                    vectorR1[k] <- 1-(igraph::compare(comrR1, comRandom1, 
                                                   method=measure))
                    vectorR2[k] <- 1-(igraph::compare(comrR2, comRandom2, 
                                                   method=measure))
                }
                
                R11[count2] <- vectorR1[k]
                R22[count2] <- vectorR2[k]
                diffR <- igraph::difference(graphRandom, Random)
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
                    measReal11[count2] <- vector1[k]
                    measReal22[count2] <- vector2[k]
                    #Random
                    graphRewireRandom <- rewireOnl(data=Random,
                                                    number=round(0.01*z))
                    comrR1 <- membershipCommunities(graph=Random, 
                                                    method=method1,
                                                    FUN=FUN1,
                                                    directed=directed,
                                                    weights=weights,
                                                    steps=steps, 
                                                    spins=spins, 
                                                    e.weights=e.weights, 
                                                    v.weights=v.weights, 
                                                    nb.trials=nb.trials)
                    comrR2 <- membershipCommunities(graph=Random, 
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
                        vectorR1[k] <- igraph::compare(comrR1, comRandom1, 
                                                       method=measure)
                        vectorR2[k] <- igraph::compare(comrR2, comRandom2,
                                                       method=measure)
                    } else {
                        vectorR1[k] <- 1-(igraph::compare(comrR1, comRandom1, 
                                                       method=measure))
                        vectorR2[k] <- 1-(igraph::compare(comrR2, comRandom2,
                                                       method=measure))
                    }
                    
                    R11[count2] <- vectorR1[k]
                    R22[count2] <- vectorR2[k]
                }
                Mean11[s] <- mean(measReal11)
                Mean22[s] <- mean(measReal22)
                MeanRandom11[s] <- mean(R11)
                MeanRandom22[s] <- mean(R22)
            }
            graph <- igraph::intersection(graph, graphRewire)
            Mean1 <-cbind(Mean1,Mean11)
            Mean2 <-cbind(Mean2,Mean22)
            graphRandom <- igraph::intersection(graphRandom, Random)
            MeanRandom1 <-cbind(MeanRandom1,MeanRandom11)
            MeanRandom2 <-cbind(MeanRandom2,MeanRandom22)
            z1 <- igraph::gsize(graph)
            print(z1)
        }
    }
    colnames(Mean1) <- nRewire 
    colnames(Mean2) <- nRewire 
    colnames(MeanRandom1) <- nRewire
    colnames(MeanRandom2) <- nRewire
    nn <- rep(nRewire, each=nrep) 
    
    ratios1 <- log2((Mean1+0.001)/(MeanRandom1+0.001))
    ratios2 <- log2((Mean2+0.001)/(MeanRandom2+0.001))
    ratios1vs2 <- log2((Mean1+0.001)/(Mean2+0.001))
    bats1 <- as.vector(ratios1)
    bats2 <- as.vector(ratios2)
    bats1vs2 <- as.vector(ratios1vs2)
    names(bats1) <- nn
    names(bats2) <- nn
    names(bats1vs2) <- nn
    res1 <- cbind(ID="ratios", t(bats1))
    res2 <- cbind(ID="ratios", t(bats2))
    res1vs2 <- cbind(ID="ratios", t(bats1vs2))
    
    output <- list(Mean1=Mean1,
                   Mean2=Mean2,
                   MeanRandom1=MeanRandom1,
                   MeanRandom2=MeanRandom2,
                   ratios1=res1,
                   ratios2=res2,
                   ratios1vs2=res1vs2)
    return(output)
}


################### PLOT COMPARISON ##################
#' plotRobinCompare
#'
#' @param graph The output of prepGraph.
#' @param model1 The Mean1 output of the comparison function.
#' @param modelR1 The MeanRandom1 output of the comparison function.
#' @param model2 The Mean2 output of the comparison function.
#' @param modelR2 The MeanRandom2 output of the comparison function.
#' @param legend The legend for the graph. The default is c("real data", 
#' "null model").
#' @param legend1vs2 The legend for the two methods.
#' The default is c("method1","method2")
#' @param title1 The title for the plot of the first method. 
#' The default is "Method1".
#' @param title2 The title for the plot of the second method.
#' The default is "Method2".
#' @param title1vs2 The title for the plot combined. The default is "Method1 vs
#' Method2".
#'
#' @return A ggplot object in a grid format.
#' @importFrom gridExtra grid.arrange
#' @export
#'
#' @examples 
#' plotRobinCompare(graph=graph, model1=Comp$Mean1, 
#' model2=Comp$Mean2,modelR1=Comp$MeanRandom1, modelR2=Comp$MeanRandom2,
#' legend=c("real data", "null model"),legend1vs2=c("Louvain", "Fast Greedy"),
#' title1="Louvain",title2="Fast Greedy",
#' title1vs2="Louvain vs Fast Greedy")
#' 
plotRobinCompare <- function(graph, model1, modelR1, model2, modelR2,
                            measure= c("vi", "nmi","split.join",
                                       "adjusted.rand"),
                            legend=c("real data", "null model"),
                            legend1vs2=c("method1", "method2"),
                            title1="Method1",title2="Method2",
                            title1vs2="Method1 vs Method2")
{
    plot1 <- plotRobin(graph=graph, model=model1, 
                    modelR=modelR1, measure=measure, legend=legend)+
                    ggtitle(title1)
    plot2 <- plotRobin(graph=graph, model=model2, modelR=modelR2, 
                       measure=measure, legend=legend)+ggtitle(title2)
    plot3 <- plotRobin(graph=graph, model=model1, modelR=model2, 
                       measure=measure,legend=legend1vs2)+ggtitle(title1vs2)
    
    gridExtra::grid.arrange(plot1,plot2,plot3, nrow=2)
}

########################### TEST ROBIN ###################

################ CREATE ITPSpline ##################
#' createITPSplineResult
#' 
#' @description creates an fdatest::ITP2 class object 
#' 
#' @param graph The output of prepGraph.
#' @param model1 The Mean output of the robinProc function (or the Mean1 
#' output of the comparison function).
#' @param model2 The MeanRandom output of the robinProc function (or the 
#' Mean2 output of the comparison function).
#' @param muParam the mu parameter for ITP2bspline (default 0).
#' @param orderParam the order parameter for ITP2bspline (default 4).
#' @param nKnots the nknots parameter for ITP2bspline (default 7).
#' @param BParam the B parameter for ITP2bspline (default 10000).
#' @param isPaired the paired parameter for ITP2bspline (default TRUE).
#'
#' @return an ITP2 object
createITPSplineResult <- function(graph, model1, model2,
                                muParam=0, orderParam=4, nKnots=7, 
                                BParam=10000, isPaired=TRUE) 
{
    nOrder <- log2(igraph::gorder(graph))
    modeled1 <- as.matrix(model1)/nOrder
    modeled2 <- as.matrix(model2)/nOrder
    ## Two populations Interval Testing Procedure with B-spline basis
    ITPresult <- fdatest::ITP2bspline(modeled1, modeled2, mu=muParam, 
                                    order=orderParam, nknots=nKnots, 
                                    B=BParam, paired=isPaired)
    return(ITPresult)
}

######################GAUSSIN PROCESS#################

#' robinGPTest
#'
#' @description This function makes a test between the curves, calculating the 
#' Bayes factor.
#' @param ratio The ratios output of the robinProc function (or the ratios1vs2 
#' output of the comparison function). 
#'
#' @return The Bayes factor
#' @import gprege
#' @export
#'
#' @examples robinGPTest(ratio=Proc$ratios)
robinGPTest <- function(ratio)
{
    gpregeOptions <- list(indexRange=(1:2), explore=FALSE, 
                         exhaustPlotRes=30, exhaustPlotLevels=10, 
                         exhaustPlotMaxWidth=100, iters=100, 
                         labels=rep(FALSE,2), display=FALSE)
    MA <- as.matrix(ratio[,c(2:dim(ratio)[2])])
    MA <- t(MA)
    vt <- unique(colnames(MA))
    ntimes <- length(vt)
    stdv <- NULL
    varv <- NULL
    for (i in c(1:ntimes)){    #ntime number of percentuage of rewire
        ind <- which(colnames(MA)==vt[i])
        stdv[i] <- sd(MA[1,ind])
        varv[i] <- var(MA[1,ind])
    }
    #sigmaest=mean(stdv)Order of the B-spline basis expansion.
    GlobalVar <- var(MA[1,])
    SigNoise <- mean(varv)/GlobalVar
    if (SigNoise>1) SigNoise <- 1
    #SigNoise=1-var(MA[2,])
    sigmaest <- 1-SigNoise
    #mod='08'
    #ntimes=50
    gpregeOptions$inithypers <- matrix( c(
        1/1000,	0,	1
        ,1/ntimes,	sigmaest, SigNoise
    ), ncol=3, byrow=TRUE)
    #gpregeOptions$inithypers <- matrix( c(
    # 1/1000,  0,	1
    #,1/ntimes,	0.8,0.2
    #), ncol=3, byrow=TRUE)
    dvet <- data.matrix(as.numeric(colnames(ratio)[-1]))
    dd <- t(data.matrix(as.numeric((ratio)[-1])))
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
#'@description The function implements the Interval Testing Procedure for 
#'testing the difference between the two curves.
#'
#' @param graph The output of prepGraph.
#' @param model1 The Mean output of the robinProc function (or the Mean1 
#' output of the comparison function).
#' @param model2 The MeanRandom output of the robinProc function (or the 
#' Mean2 output of the comparison function).
#' @param legend The legend for the graph. The default is c("real data", 
#' "null model").
#'
#' @return Two plots.
#' @import igraph ggplot2 
#' @importFrom fdatest ITP2bspline
#' @export
#'
#' @examples robinFDATest(graph=graph, model1=Proc$Mean,
#' model2=Proc$MeanRandom)
robinFDATest <- function(graph,model1,model2, 
                        legend=c("real data", "null model"))
{
     N <- igraph::vcount(graph)
    mvimodel1 <- cbind(as.vector((model1)/log2(N)), legend[1], 
                        seq(1,10), rep((seq(0,60,5)/100), each=10))
    mvimodel2 <- cbind(as.vector((model2)/log2(N)), legend[2], 
                        seq(11,20),rep((seq(0,60,5)/100),each=10))
    mvi <- rbind(mvimodel1, mvimodel2)
    colnames(mvi) <- c("mvi","model","s","percPert")
    dataFrame <- data.frame(mvi)
    plot1 <- ggplot2::ggplot(dataFrame, ggplot2::aes(x=percPert, 
                y=as.numeric(as.character(mvi)), color= model, group=s)) + 
        ggplot2::geom_line() + 
        ggplot2::xlab("Percentage of perturbation") +
        ggplot2::ylab("Measure")+
        ggplot2::ggtitle("Robin plot")
    plot1
    print(plot1)
   
    perc <- rep((seq(0,60,5)/100))
    ITPresult <- createITPSplineResult(graph, model1, model2)
    plot2 <- plot(ITPresult, main='Measure', xrange=c(0,0.6), xlab='perturbation', 
                ylab="Measure")
    lines(perc, rep(0.05, 13), type="l", col="red")
    
    print(plot2)
}  

########### AREA UNDER THE CURVE    ##############
#' robinAUCTest
#'
#' @description Calculate the area under the two curves with a spline approach. 
#' @param graph The output of prepGraph.
#' @param model1 The Mean output of the robinProc function (or the Mean1 
#' output of the comparison function).
#' @param model2 The MeanRandom output of the robinProc function (or the 
#' Mean2 output of the comparison function).
#'
#' @return A list
#' @importFrom DescTools AUC
#' @importFrom igraph vcount
#' @export
#'
#' @examples 
#' robinAUCTest(graph=graph, model1=Proc$Mean, model2=Proc$MeanRandom)
robinAUCTest <- function(graph,model1,model2)
{
    N <- igraph::vcount(graph)
    mvimeanmodel1 <- cbind(as.vector((apply(model1, 2, mean))/log2(N)))
    mvimeanmodel2 <- cbind(as.vector((apply(model2, 2, mean))/log2(N)))
    area1 <- DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel1, 
                          method ="spline")
    area2 <- DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel2, 
                          method ="spline")
    output <- list(area1=area1,area2=area2)
return(output)
}

########  ALL CURVES TOGETHER ######
# allCurves <- function(graph,graphRandom,
#                       FUN1=NULL,
#                       FUN2=NULL,
#                       type=c("dependent", "independent"),
#                       directed=FALSE,
#                       weights=NULL, 
#                       steps=4, 
#                       spins=25, 
#                       e.weights=NULL, 
#                       v.weights=NULL, 
#                       nb.trials=10)
# {
#     seq <- c("walktrap", "edgeBetweenness", 
#              "fastGreedy", "louvain", "spinglass", 
#              "leadingEigen", "labelProp", "infomap",
#              "optimal", "other")
#     for(i in 1:9) {
#         comp <- comparison(graph,graphRandom, 
#                                   method1=seq[i],
#                                   method2=seq[i+1],
#                                   FUN1=NULL,
#                                   FUN2=NULL,
#                                   type=c("dependent", "independent"),
#                                   directed=FALSE,
#                                   weights=NULL, 
#                                   steps=4, 
#                                   spins=25, 
#                                   e.weights=NULL, 
#                                   v.weights=NULL, 
#                                   nb.trials=10)
#         plotRobinCompare(graph=graph, model1=comp$viMean1, model2=comp$viMean2,
#                          modelR1=comp$viMeanRandom1, modelR2=comp$viMeanRandom2,
#                          legend=c("real data", "null model"),
#                          legend1vs2=c(seq[i],seq[i+1]),
#                          title1=seq[i], title2=seq[i+1],
#                          title1vs2=c(seq[i],seq[i+1]))
#                          }
#         
# }   


