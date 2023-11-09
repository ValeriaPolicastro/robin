####### GRAPH RANDOM WEIGHTED#########
#' randomWeight
#'
#' @description This function randomly rewires the edges while preserving the original graph's
#' degree distribution.
#' @param graph The output of prepGraph.
#' @param distrib Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "NegBinom" for negative binomial.
#' @param verbose flag for verbose output (default as FALSE)
#'
#' @return An igraph object, a randomly rewired graph.
#' @import igraph
#' @export
#'
#' @examples
#' my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' E(graph)$weight <- round(runif(ecount(graph),min=1,max=10))
#' graphRandom <- randomWeight(graph=graph)
randomWeight <- function(graph, distrib = "NegBinom",verbose=FALSE)
{
    if(verbose) cat("Randomizing the graph edges.\n")
    v <- igraph::vcount(graph) ## number of vertex
    numberPerturbAll <- round((v*(v-1))/2, 0)
    adj <- as_adjacency_matrix(graph, attr="weight", sparse = FALSE)
    graphRandom <- as.matrix(perturbR::rewireR(sym.matrix=adj, 
                                               nperturb=numberPerturbAll, 
                                               dist = distrib))
    #rewiring for z all the edges
    graphRandom <- graph_from_adjacency_matrix(graphRandom,weighted = TRUE, 
                                               mode="undirected")
    return(graphRandom)
}






########## ROBIN COMPARE WEIGHTED#############

#' robinCompareFastWeight
#'
#' @description This function compares two community detection algorithms, from
#' weighted networks. Is the parallelized and faster version.
#' @param graph The output of prepGraph.
#' @param method1 The first clustering method, one of "walktrap",
#' "edgeBetweenness", "fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal","leiden".
#' @param method2 The second custering method one of "walktrap",
#' "edgeBetweenness","fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap","optimal","leiden".
#' @param measure The stability measure, one of "vi", "nmi", "split.join",
#' "adjusted.rand" all normalized and used as distances.
#' "nmi" refers to 1- nmi and "adjusted.ran" refers to 1-adjusted.rand.
#' @param ncores number of CPU cores to use.(default is 2) For a faster
#' execution we suggest to use ncores=(detectCores(logical = FALSE)-1)
#' @param FUN1 personal designed function when method1 is "others".
#' see \code{\link{methodCommunity}}.
#' @param FUN2 personal designed function when method2 is "others".
#' see \code{\link{methodCommunity}}.
#' @param weights This argument is not settable for "infomap" method.
#' @param steps This argument is settable only for "leadingEigen"and"walktrap"
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method.
#' @param objective_function Whether to use the Constant Potts Model (CPM) or 
#' modularity. Must be either "CPM" or "modularity".
#' @param n_iterations the number of iterations to iterate the Leiden algorithm. 
#' Each iteration may improve the partition further.This argument is settable 
#' only for "leiden".
#' @param resolution only for "louvain" and "leiden". Optional resolution 
#' parameter, lower values typically yield fewer, larger clusters (default=1).
#' @param verbose flag for verbose output (default as TRUE).
#' @param distrib Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "NegBinom" for negative binomial.
#'
#' @return A list object with two matrices:
#' - the matrix "Mean1" with the means of the procedure for the first method
#' - the matrix "Mean2" with the means of the procedure for the second method
#'
#' @import igraph parallel perturbR
#' @export
#'
#' @examples my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' E(graph)$weight <- round(runif(ecount(graph),min=1,max=10))
#' robinCompareFastWeight(graph=graph, method1="louvain",
#' method2="fastGreedy", measure="vi")


robinCompareFastWeight <- function(graph,
                                   method1=c("walktrap", "edgeBetweenness", "fastGreedy",
                                             "leadingEigen","louvain","spinglass",
                                             "labelProp","infomap","optimal","leiden",
                                             "other"),
                                   method2=c("walktrap", "edgeBetweenness", "fastGreedy",
                                             "leadingEigen","louvain","spinglass",
                                             "labelProp","infomap","optimal","leiden",
                                             "other"),
                                   measure= c("vi", "nmi","split.join", "adjusted.rand"),
                                   ncores=2,
                                   FUN1=NULL, FUN2=NULL,
                                   directed=FALSE, weights=NULL, steps=4,
                                   spins=25, e.weights=NULL, v.weights=NULL,
                                   nb.trials=10,resolution = 1,n_iterations=2,
                                   objective_function = c("CPM", "modularity"),
                                   verbose=TRUE, distrib="NegBinom")
{
    method1 <- match.arg(method1)
    method2 <- match.arg(method2)
    measure <- match.arg(measure)
    comReal1 <- membershipCommunities(graph=graph, method=method1,
                                      FUN=FUN1,
                                      directed=directed,
                                      weights=weights,
                                      steps=steps,
                                      spins=spins,
                                      e.weights=e.weights,
                                      v.weights=v.weights,
                                      nb.trials=nb.trials,
                                      resolution=resolution,
                                      objective_function = objective_function,
                                      n_iterations=n_iterations)
    comReal2 <- membershipCommunities(graph=graph, method=method2,
                                      FUN=FUN2,
                                      directed=directed,
                                      weights=weights,
                                      steps=steps,
                                      spins=spins,
                                      e.weights=e.weights,
                                      v.weights=v.weights,
                                      nb.trials=nb.trials,
                                      resolution=resolution,
                                      objective_function = objective_function,
                                      n_iterations=n_iterations)
    
    N <- igraph::vcount(graph)
    de <- round((N*(N-1))/2, 0)
    Measure <- NULL
    vector1 <- NULL
    vector2 <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detecting robin method independent type, wait it can take time it depends on the size of the network.\n")
    vet1 <- seq(5, 60, 5)
    vet <- round(vet1*de/100, 0)
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl,varlist =c("graph","method1","method2","directed",
                                          "weights","steps","spins", "e.weights",
                                          "v.weights", "nb.trials","measure","comReal1",
                                          "comReal2","N","verbose","FUN1","FUN2","distrib",
                                          "resolution","objective_function",
                                          "n_iterations"),
                            envir=environment())
    zlist <- parallel::clusterApply(cl,vet, function(z)
    {
        
        
        
        MeansList <- lapply(1:5, function(s)
        {
            
            
            adj <- igraph::as_adjacency_matrix(graph, attr="weight", sparse = FALSE)
            gR <- as.matrix(perturbR::rewireR(adj, z,dist = distrib))
            graphRList <- igraph::graph_from_adjacency_matrix(gR,weighted = TRUE,
                                                              mode="undirected")
            
            
            comr1 <- robin::membershipCommunities(graph=graphRList,
                                                  method=method1,
                                                  FUN=FUN1,
                                                  directed=directed,
                                                  weights=weights,
                                                  steps=steps,
                                                  spins=spins,
                                                  e.weights=e.weights,
                                                  v.weights=v.weights,
                                                  nb.trials=nb.trials,
                                                  resolution=resolution,
                                                  objective_function = objective_function,
                                                  n_iterations=n_iterations)
            
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
                                                  nb.trials=nb.trials,
                                                  resolution=resolution,
                                                  objective_function = objective_function,
                                                  n_iterations=n_iterations)
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
    parallel::stopCluster(cl)
    Measure1 <- do.call(cbind, lapply(zlist, function(z) z$Measure1))
    Measure2 <- do.call(cbind, lapply(zlist, function(z) z$Measure2))
    
    Measure1 <- cbind(rep(0, 5), Measure1)
    Measure2 <- cbind(rep(0, 5), Measure2)
    
    colnames(Measure1) <- nRewire
    colnames(Measure2) <- nRewire
    return(list(Mean1=Measure1,
                Mean2=Measure2))
}



#' rewireComplWeight
#'
#' @description rewires the weighted graph, creates the communities and
#' compares the communities through different measures.
#'
#' @param data The output of prepGraph.
#' @param number Number of rewiring trials to perform.
#' @param community Community to compare with.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness",
#' "fastGreedy", "louvain", "spinglass","leadingEigen", "labelProp", "infomap",
#' "optimal","leiden", "other".
#' @param FUN see \code{\link{methodCommunity}}.
#' @param measure The measure for the comparison of the communities "vi", "nmi",
#' "split.join", "adjusted.rand".
#' @param weights this argument is not settable for "infomap" method.
#' @param steps this argument is settable only for "leadingEigen"and"walktrap"
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method
#' @param objective_function Whether to use the Constant Potts Model (CPM) or 
#' modularity. Must be either "CPM" or "modularity".
#' @param n_iterations the number of iterations to iterate the Leiden algorithm. 
#' Each iteration may improve the partition further.This argument is settable 
#' only for "leiden".
#' @param resolution only for "louvain" and "leiden". Optional resolution
#'  parameter, lower values typically yield fewer, larger clusters (default=1).
#' @param distrib Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "NegBinom" for negative binomial.
#' @keywords internal
#'
rewireComplWeight <- function(data, number, community,
                              method=c("walktrap", "edgeBetweenness",
                                       "fastGreedy", "louvain", "spinglass",
                                       "leadingEigen", "labelProp", "infomap",
                                       "optimal","leiden", "other"),
                              FUN=NULL,
                              measure= c("vi", "nmi","split.join", "adjusted.rand"),
                              directed=FALSE, weights=NULL, steps=4, spins=25,
                              e.weights=NULL, v.weights=NULL, nb.trials=10, 
                              resolution = 1,n_iterations=2,
                              objective_function = c("CPM", "modularity"),
                              distrib="NegBinom")
{
    method <- match.arg(method)
    measure <- match.arg (measure)
    adj <- as_adjacency_matrix(data, attr="weight", sparse = FALSE)
    graphRewire <- as.matrix(perturbR::rewireR(adj, number,dist = distrib))
    graphRewire <- graph_from_adjacency_matrix(graphRewire,weighted = TRUE, mode="undirected")
    comR <- membershipCommunities(graph=graphRewire, method=method, FUN=FUN,
                                  directed=directed, weights=weights, steps=steps,
                                  spins=spins, e.weights=e.weights,
                                  v.weights=v.weights, nb.trials=nb.trials,
                                  resolution=resolution,
                                  objective_function = objective_function,
                                  n_iterations=n_iterations)
    Measure <- igraph::compare(community, comR, method=measure)
    output <- list(Measure=Measure, graphRewire=graphRewire)
    
    return(output)
}

########### ROBIN ROBUST WEIGHTED #########

#' robinRobustWeighted
#'
#' @description This functions implements a procedure to examine the stability
#' of the partition recovered by some algorithm against random perturbations
#' of the original graph structure for weighted network.
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
#' @param weights this argument is not settable for "infomap" method.
#' @param steps this argument is settable only for "leadingEigen"and"walktrap"
#' method.
#' @param spins This argument is settable only for "infomap" method.
#' @param e.weights This argument is settable only for "infomap" method.
#' @param v.weights This argument is settable only for "infomap" method.
#' @param nb.trials This argument is settable only for "infomap" method.
#' @param directed This argument is settable only for "edgeBetweenness" method.
#' @param objective_function Whether to use the Constant Potts Model (CPM) or 
#' modularity. Must be either "CPM" or "modularity".
#' @param n_iterations the number of iterations to iterate the Leiden algorithm. 
#' Each iteration may improve the partition further.This argument is settable 
#' only for "leiden".
#' @param resolution only for "louvain" and "leiden". Optional resolution 
#' parameter, lower values typically yield fewer, larger clusters (default=1).
#' @param distrib Option to rewire in a manner that retains overall graph weight 
#' regardless of distribution of edge weights. This option is invoked by putting 
#' any text into this field. Defaults to "NegBinom" for negative binomial.
#' @param verbose flag for verbose output (default as TRUE).
#'
#' @return A list object with two matrices:
#' - the matrix "Mean" with the means of the procedure for the graph
#' - the matrix "MeanRandom" with the means of the procedure for the random graph.
#'
#' @import igraph perturbR
#' @export
#' @examples my_file <- system.file("example/football.gml", package="robin")
#' graph <- prepGraph(file=my_file, file.format="gml")
#' E(graph)$weight <- round(runif(ecount(graph),min=1,max=10))
#' graphRandom <- randomWeight(graph=graph)
#' robinRobustWeighted(graph=graph, graphRandom=graphRandom,
#' method="louvain", measure="vi")
robinRobustWeighted <- function(graph, graphRandom,
                                method=c("walktrap", "edgeBetweenness",
                                         "fastGreedy", "louvain", "spinglass",
                                         "leadingEigen", "labelProp", "infomap",
                                         "optimal","leiden", "other"),
                                FUN=NULL, measure= c("vi", "nmi","split.join", 
                                                     "adjusted.rand"),
                                directed=FALSE, weights=NULL,
                                steps=4, spins=25, e.weights=NULL, v.weights=NULL,
                                nb.trials=10, resolution = 1,n_iterations=2,
                                objective_function = c("CPM", "modularity"),
                                distrib="NegBinom",
                                verbose=TRUE)
{
    measure <- match.arg(measure)
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
                                     nb.trials=nb.trials,
                                     resolution=resolution,
                                     objective_function = objective_function,
                                     n_iterations=n_iterations) # real network
    
    comRandom <- membershipCommunities(graph=graphRandom, method=method,
                                       FUN=FUN,
                                       directed=directed,
                                       weights=weights,
                                       steps=steps,
                                       spins=spins,
                                       e.weights=e.weights,
                                       v.weights=v.weights,
                                       nb.trials=nb.trials,
                                       resolution=resolution,
                                       objective_function = objective_function,
                                       n_iterations=n_iterations) # random network
    #stopifnot(length(table(comRandom))>1)
    #if(length(table(comRandom))==1) {stop("Not random graph")}
    N <- igraph::vcount(graph)
    de <- round((N*(N-1))/2, 0)
    Measure <- NULL
    vector <- NULL
    vectRandom <- NULL
    graphRewireRandom <- NULL
    graphRewire <- NULL
    count <- 1
    nRewire <- seq(0,60,5)
    if(verbose) cat("Detected robin method \n")
    
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
            Real <- rewireComplWeight(data=graph,
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
                                      nb.trials=nb.trials,
                                      resolution=resolution,
                                      objective_function = objective_function,
                                      n_iterations=n_iterations,
                                      distrib=distrib)
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
            Random <- rewireComplWeight(data=graphRandom,
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
                                        nb.trials=nb.trials,
                                        resolution=resolution,
                                        objective_function = objective_function,
                                        n_iterations=n_iterations,
                                        distrib=distrib)
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
                Real <- rewireComplWeight(data=graphRewire,
                                          number=round(0.01*de),
                                          community=comReal,
                                          method=method,
                                          measure=measure,
                                          directed=directed,
                                          weights=weights,
                                          steps=steps,
                                          spins=spins,
                                          e.weights=e.weights,
                                          v.weights=v.weights,
                                          nb.trials=nb.trials,
                                          resolution=resolution,
                                          objective_function = objective_function,
                                          n_iterations=n_iterations,
                                          distrib=distrib)
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
                Random <- rewireComplWeight(data=graphRewireRandom,
                                            number=round(0.01*de),
                                            community=comRandom,
                                            method=method,
                                            measure=measure,
                                            directed=directed,
                                            weights=weights,
                                            steps=steps,
                                            spins=spins,
                                            e.weights=e.weights,
                                            v.weights=v.weights,
                                            nb.trials=nb.trials,
                                            resolution=resolution,
                                            objective_function = objective_function,
                                            n_iterations=n_iterations,
                                            distrib=distrib)
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

########## READ NSL ###########
read_nsl=function(file){
    df=read.table(file,skip=1)
    colnames(df)=c('V1','V2','weight')
    gg=graph_from_data_frame(df,directed = FALSE)
    # gg <- set_edge_attr(gg, "weight", value= df$weights)
    return(gg)
}
