
######PREPARATION GRAPH########## 

#' prepGraph
#' 
#' @description The prepGraph function is able to read graphs from a file and 
#' to prepare them for the analysis.
#' @param file The file to read from.
#' @param file.format Character constant giving the file format. Right now
#' as_edgelist, pajek, graphml, gml, ncol, lgl, dimacs and graphdb are 
#' supported.
#'
#' @return A simple graph
#' @import igraph
#' @export
#'
#' @examples graph <- prepGraph(file=my_file, file.format="edgelist")
prepGraph <- function(file,
                    file.format=c("edgelist", "pajek", "ncol", "lgl", "graphml",
                                    "dimacs", "graphdb", "gml", "dl"))
{
    net <- igraph::read_graph(file=file, format=file.format, directed=FALSE)
    ind <- igraph::V(net)[degree(net) == 0] #isolate node
    graph <- igraph::delete.vertices(net, ind)
    graph <- igraph::simplify(graph)
    
    ##2
    # edge <- read.table(file,colClasses = "character", quote="\"", header=header)
    # edge <- as.matrix(edge)
    # graph<-graph_from_edgelist(edge,direct=direct)
    # graph <- igraph::simplify(graph)
    
    ##1
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
    #graphs which do not contain loop and multiple edges.
    return(graph)
}


######GRAPH RANDOM#########
#' random
#'
#' @description Randomly rewire the edges while preserving the original graph's 
#' degree distribution.
#' @param graph The output of prepGraph
#'
#' @return A randomly rewired graph
#' @import igraph
#' @export
#'
#' @examples graphRandom <- random(graph=graph)
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
#' @description This function gives the membership vector of the community 
#' structure. The community structure was found by functions 
#' implemented in igraph.
#' @param graph The input graph created with prepGraph
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal"
#' @param weights this argument is not settable for "infomap" method
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" 
#' method
#' @param spins This argument is settable only for "infomap" method
#' @param e.weights This argument is settable only for "infomap" method
#' @param v.weights This argument is settable only for "infomap" method
#' @param nb.trials This argument is settable only for "infomap" method
#' @param directed This argument is settable only for "edgeBetweenness" method
#'
#' @return Membership vector of the community structure
#' @import igraph
#' @export
#'
#' @examples methodCommunity (graph=graph, method="louvain")
#' #with all the method implemented in igraph
methodCommunity <- function(graph, 
                            method,
                            directed=FALSE,
                            weights=NULL, 
                            steps=4, 
                            spins=25, 
                            e.weights=NULL, 
                            v.weights=NULL, 
                            nb.trials=10)
{
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
                                v.weights=v.weights, nb.trials=nb.trials)
    )
    return(membership(communities))
}

################ PLOT GRAPH ###############
#' plotGraph
#'
#' @description Graphical interactive representation of the network
#' @param graph The output of prepGraph
#'
#' @return An interactive plot
#' @import networkD3
#' @export
#'
#' @examples plotGraph (graph)
plotGraph <- function(graph)
{
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph)
    plot <- networkD3::simpleNetwork(graph_d3$links)
    return(plot)
}   


######################## PLOT COMMUNITIES ##############
#' plotCommu
#'
#' @description An interactive 3D plot of the communities
#' @param graph The output of prepGraph
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap",
#' "optimal"
#'
#' @return An interactive plot
#' @import networkD3
#' @export
#'
#' @examples plotCommu(graph=graph, method="louvain")
plotCommu <- function(graph, method)
{
    members<-methodCommunity(graph=graph,method=method)
    # Convert to object suitable for networkD3
    graph_d3 <- networkD3::igraph_to_networkD3(g=graph, group = members)
    # Create force directed network plot
    plot <- networkD3 ::forceNetwork(Links = graph_d3$links, Nodes = graph_d3$nodes,
                         Source ='source', 
                 Target ='target', NodeID ='name',Group ='group')
    return(plot)
}


######### REWIRE COMPLETE ########
#' rewireCompl
#' 
#' @param data The output of prepGraph
#' @param number Number of rewiring trials to perform.
#' @param community Community to compare with.
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap"
#' @param weights this argument is not settable for "infomap" method
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" 
#' method
#' @param spins This argument is settable only for "infomap" method
#' @param e.weights This argument is settable only for "infomap" method
#' @param v.weights This argument is settable only for "infomap" method
#' @param nb.trials This argument is settable only for "infomap" method
#' @param directed This argument is settable only for "edgeBetweenness" method
#'
#' @return A list object

rewireCompl <- function(data, number, community, method,
                        directed=FALSE,
                        weights=NULL, 
                        steps=4, 
                        spins=25, 
                        e.weights=NULL, 
                        v.weights=NULL, 
                        nb.trials=10)
{
    graphRewire <- igraph::rewire(data,
                                  with=keeping_degseq(loops=FALSE,niter=number))
    comR  <- methodCommunity(graph=graphRewire, method=method,
                            directed=directed,
                            weights=weights,
                            steps=steps, 
                            spins=spins, 
                            e.weights=e.weights, 
                            v.weights=v.weights, 
                            nb.trials=nb.trials)
    VI <- compare(community, comR, method="vi")
    output <- list(VI=VI, graphRewire=graphRewire)
    return(output)#non serve per tutti
}
#perturba il grafo e ricalcola le community
#compare the community through VI

#########REWIRE ONLY ###########
#' rewireOnl
#'
#' @param data The output of prepGraph
#' @param number Number of rewiring trials to perform.
#'
#' @return A graph object
rewireOnl <- function(data, number)
{
    graphRewire <- igraph::rewire(data,
                                  with=keeping_degseq(loops=FALSE,niter=number))
    return(graphRewire)
}

########  ROBIN PROCEDURE #######
#' robinProc
#'
#' @description A procedure to examine the stability of the partition recovered 
#' against random perturbations of the original graph structure.
#' @param graph The output of prepGraph
#' @param graphRandom The output of random command
#' @param type The type of robin costruction dependent or independent data
#' @param method The clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap"
#' @param weights this argument is not settable for "infomap" method
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" 
#' method
#' @param spins This argument is settable only for "infomap" method
#' @param e.weights This argument is settable only for "infomap" method
#' @param v.weights This argument is settable only for "infomap" method
#' @param nb.trials This argument is settable only for "infomap" method
#' @param directed This argument is settable only for "edgeBetweenness" method
#'
#' @return A list object
#' @import igraph
#' @export
#'
#' @examples Proc<- robinProc(graph=graph,graphRandom=graphRandom, method="louvain",type="independent")
robinProc <- function(graph, graphRandom, method,
                type=c("dependent", "independent"),
                directed=FALSE,
                weights=NULL, 
                steps=4, 
                spins=25, 
                e.weights=NULL, 
                v.weights=NULL, 
                nb.trials=10) 
{
    type <- match.arg(type)
    nrep <- 10
    comReal <- methodCommunity(graph=graph, method=method, directed=directed,
                                weights=weights, 
                                steps=steps, 
                                spins=spins, 
                                e.weights=e.weights, 
                                v.weights=v.weights, 
                                nb.trials=nb.trials) # real network

    comRandom <- methodCommunity(graph=graphRandom, method=method, 
                                directed=directed,
                                weights=weights,
                                steps=steps, 
                                spins=spins, 
                                e.weights=e.weights, 
                                v.weights=v.weights, 
                                nb.trials=nb.trials) # random network
    de <-igraph::gsize(graph)
    VI <- NULL
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
        vi <- matrix(0, nrep^2, length(nRewire))
        viRandom <- matrix(0, nrep^2, length(nRewire))
        viMeanRandom <- matrix(0, nrep, length(nRewire))
        viMean <- matrix(0, nrep, length(nRewire))
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
                vector[k] <- Real$VI
                vi[count2, count] <- vector[k]
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
                vectRandom[k] <- Random$VI
                viRandom[count2, count] <- vectRandom[k]
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
                    vector[k] <- Real$VI
                    vi[count2, count] <- vector[k]
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
                    vectRandom[k] <- Random$VI
                    viRandom[count2, count] <- vectRandom[k]
                }
                viMeanRandom[s, count] <- mean(vectRandom)
                viMean[s, count] <- mean(vector)
                print(count)
            }
        }
  #DEPENDENT 
    }else{
        z <- round((5*de)/100, 0) #the 5% of the edges
        vi <- rep(0, nrep^2)
        vi1 <- NULL
        viRandom <- rep(0, nrep^2)
        viRandom1 <- NULL
        viMeanRandom <- rep(0, nrep)
        viMean <- rep(0, nrep)
        viMeanRandom1 <- NULL
        viMean1 <- NULL
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
                comr <- methodCommunity(graph=graphRewire,
                                        method=method,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials)
                vector[k] <- compare(comReal, comr, method="vi")
                vi1[count2] <- vector[k]
                diff <-igraph::difference(graph, graphRewire)
                
                ###RANDOM
                graphRewireRandom <- rewireOnl(data=graphRandom, number=z)
                graphRewireRandom <- igraph::union(graphRewireRandom, diffR)
                comr <- methodCommunity(graph=graphRewireRandom,
                                        method=method,
                                        directed=directed,
                                        weights=weights,
                                        steps=steps, 
                                        spins=spins, 
                                        e.weights=e.weights, 
                                        v.weights=v.weights, 
                                        nb.trials=nb.trials)
                vectRandom[k] <- compare(comRandom, comr, method="vi")
                viRandom1[count2] <- vectRandom[k]
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
                    vector[k] <- Real$VI
                    vi1[count2] <- vector[k]
                    
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
                    vectRandom[k] <- Random$VI
                    viRandom1[count2] <- vectRandom[k]
                }
                viMean1[s] <- mean(vi1)
                viMeanRandom1[s] <- mean(viRandom1)  
            }
            graph <- igraph::intersection(graph, graphRewire)
            graphRandom <- igraph::intersection(graphRandom, graphRewireRandom)
            viRandom <- cbind(viRandom,viRandom1)
            vi <- cbind(vi,vi1)
            viMean <- cbind(viMean,viMean1)
            viMeanRandom <- cbind(viMeanRandom,viMeanRandom1)
            z1 <- igraph::gsize(graph)
            print(z1)
            #z2 <- igraph::gsize(graphRandom)
        }
    }
    colnames(viRandom) <- nRewire
    colnames(vi) <- nRewire
    colnames(viMeanRandom) <- nRewire
    colnames(viMean) <- nRewire
    nn <- rep(nRewire, each=nrep) 
    ratios <- log2((viMean+0.001)/(viMeanRandom+0.001))
    #rapporto tra la media delle distanze VI tra il modello reale e quello
    ##perturbato e la media delle distanze tra il random e la sua perturbazione
    bats <- as.vector(ratios)
    names(bats) <- nn
    resBats <- cbind(ID="ratios", t(bats))#la trasposta del rapporto
    
    output <- list(vi=vi,
                    viRandom=viRandom,
                    viMean=viMean,
                    viMeanRandom=viMeanRandom,
                    ratios=resBats
                    )
      return(output)
    
    
}



############PLOT##############
#' plotRobin
#'
#' @description The plot of the two VI curves, the null model and 
#' the real graph.
#' @param graph The output of prepGraph
#' @param model The viMean output of the robinProc function or the viMean output 
#' of the comparison function. The default value is viMean.
#' @param modelR The viMeanRandom output of the iter function or the viMean2
#'  output of the comparison function. The default value is viMean.
#' @param legend The legend for the graph. The default is c("real data", 
#' "null model")
#'
#' @return A plot
#' @export
#'
#' @examples 
plotRobin <- function(graph,
                      model=Proc$viMean,
                      modelR=Proc$viMeanRandom, 
                      legend=c("real data", "null model"))
{   
    N <- igraph::vcount(graph)
    mvimodel1 <- cbind(as.vector((apply(model, 2, mean))/log2(N)),legend[1])
    mvimodel2 <- cbind(as.vector((apply(modelR, 2, mean))/log2(N)),legend[2])
    percPert <- rep((seq(0,60,5)/100),2)
    mvi <- rbind(mvimodel1,mvimodel2)
    colnames(mvi) <- c("mvi","model")
    dataFrame <- data.frame(percPert,mvi)
    plot <- ggplot2::  ggplot(dataFrame, aes(x=percPert, 
                                             y=as.numeric(as.character(mvi)), 
                                             colour=model,group=factor(model)))+ 
        geom_line()+
        geom_point()+ 
        xlab("Percentage of perturbation") +
        ylab("Variation of Information (VI)")+
        ggtitle("Robin plot")
    plot+ scale_y_continuous(limits=c(0,0.6))
}

############### COMPARISON DIFFERENT METHODS ##########
####con la stessa perturbazione vengono testati due metodi
#' comparison
#'
#' @param graph The input graph prepNet
#' @param method1 The first clustering method, one of "walktrap", 
#' "edgeBetweenness", "fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap"
#' @param method2 The second custering method one of "walktrap",
#' "edgeBetweenness","fastGreedy", "louvain", "spinglass", "leadingEigen",
#' "labelProp", "infomap"
#' @param type The type of robin costruction dependent or independent data
#' @param weights this argument is not settable for "infomap" method
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" 
#' method
#' @param spins This argument is settable only for "infomap" method
#' @param e.weights This argument is settable only for "infomap" method
#' @param v.weights This argument is settable only for "infomap" method
#' @param nb.trials This argument is settable only for "infomap" method
#' @param directed This argument is settable only for "edgeBetweenness" method
#' @param graphRandom The randomly rewired graph
#'
#' @return A list object
#' @export
#'
#' @examples

comparison <- function(graph,graphRandom, 
                       method1=c("walktrap", "edgeBetweenness", "fastGreedy",
                                 "leadingEigen","louvain","spinglass",
                                 "labelProp","infomap"),
                       method2=c("walktrap", "edgeBetweenness", "fastGreedy",
                                 "leadingEigen","louvain","spinglass",
                                 "labelProp","infomap"),
                       type=c("dependent", "independent"),
                       directed=FALSE,
                       weights=NULL, 
                       steps=4, 
                       spins=25, 
                       e.weights=NULL, 
                       v.weights=NULL, 
                       nb.trials=10)
    {
    type <- match.arg(type)
    nrep <- 10
    comReal1 <- methodCommunity(graph=graph, method=method1,directed=directed,
                               weights=weights,
                               steps=steps, 
                               spins=spins, 
                               e.weights=e.weights, 
                               v.weights=v.weights, 
                               nb.trials=nb.trials) 
    comReal2 <- methodCommunity(graph=graph, method=method2,directed=directed,
                                weights=weights,
                                steps=steps, 
                                spins=spins, 
                                e.weights=e.weights, 
                                v.weights=v.weights, 
                                nb.trials=nb.trials)
    comRandom1 <- methodCommunity(graph=graphRandom, method=method1, 
                                 directed=directed,
                                 weights=weights,
                                 steps=steps, 
                                 spins=spins, 
                                 e.weights=e.weights, 
                                 v.weights=v.weights, 
                                 nb.trials=nb.trials) 
    comRandom2 <- methodCommunity(graph=graphRandom, method=method2, 
                                 directed=directed,
                                 weights=weights,
                                 steps=steps, 
                                 spins=spins, 
                                 e.weights=e.weights, 
                                 v.weights=v.weights, 
                                 nb.trials=nb.trials) 
    de <- igraph::gsize(graph)
    VI <- NULL
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
        vi1 <- matrix(0, nrep^2, length(nRewire))
        vi2 <- matrix(0, nrep^2, length(nRewire))
        viR1 <- matrix(0, nrep^2, length(nRewire))
        viR2 <- matrix(0, nrep^2, length(nRewire))
        viMean1 <- matrix(0, nrep, length(nRewire))
        viMean2 <- matrix(0, nrep, length(nRewire))
        viMeanRandom1 <- matrix(0, nrep, length(nRewire))
        viMeanRandom2 <- matrix(0, nrep, length(nRewire))
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
                comr1 <- methodCommunity(graph=graphRewire, method=method1,
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                comr2 <- methodCommunity(graph=graphRewire, method=method2, 
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                vector1[k] <- compare(comr1, comReal1, method="vi")
                vector2[k] <- compare(comr2, comReal2, method="vi")
                vi1[count2, count] <- vector1[k]
                vi2[count2, count] <- vector2[k]
                Random <- rewireOnl(data=graphRandom, number=z)
                comrR1 <- methodCommunity(graph=Random, method=method1,
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                comrR2 <- methodCommunity(graph=Random, method=method2, 
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                vectorR1[k] <- compare(comrR1, comRandom1, method="vi")
                vectorR2[k] <- compare(comrR2, comRandom2, method="vi")
                viR1[count2, count] <- vectorR1[k]
                viR2[count2, count] <- vectorR2[k]
                for(k in c(2:nrep))
                {
                    count2 <- count2+1
                    graphRewire <- rewireOnl(data=graphRewire,
                                            number=round(0.01*z))
                    comr1 <- methodCommunity(graph=graphRewire, method=method1,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    comr2 <- methodCommunity(graph=graphRewire, method=method2,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    vector1[k] <- compare(comr1, comReal1, method="vi")
                    vector2[k] <- compare(comr2, comReal2, method="vi")
                    vi1[count2, count] <- vector1[k]
                    vi2[count2, count] <- vector2[k]
                    
                    graphRewireRandom <- rewireOnl(data=Random,
                                             number=round(0.01*z))
                    comrR1 <- methodCommunity(graph=Random, method=method1,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    comrR2 <- methodCommunity(graph=Random, method=method2,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    vectorR1[k] <- compare(comrR1, comRandom1, method="vi")
                    vectorR2[k] <- compare(comrR2, comRandom2, method="vi")
                    viR1[count2, count] <- vectorR1[k]
                    viR2[count2, count] <- vectorR2[k]
                   
                }
                viMean1[s, count] <- mean(vector1)
                viMean2[s, count] <- mean(vector2)
                viMeanRandom1[s, count] <- mean(vectorR1)
                viMeanRandom2[s, count] <- mean(vectorR2)
            }
        }
        
}else{
        z <- round((5*de)/100, 0)
        vi1 <- rep(0, nrep^2)
        vi11 <- NULL
        viR11<-NULL
        vi2 <- rep(0, nrep^2)
        vi22 <- NULL
        viR22 <- NULL
        viMean1 <- rep(0, nrep)
        viMean2 <-rep(0, nrep)
        viMean11 <- NULL
        viMean22 <- NULL
        viMeanRandom1 <- rep(0, nrep)
        viMeanRandom2 <-rep(0, nrep)
        viMeanRandom11 <- NULL
        viMeanRandom22 <- NULL
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
                comr1 <- methodCommunity(graph=graphRewire, method=method1,
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                comr2 <- methodCommunity(graph=graphRewire, method=method2,
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                vector1[k] <- compare(comr1, comReal1, method="vi")
                vector2[k] <- compare(comr2, comReal2, method="vi")
                vi11[count2] <- vector1[k]
                vi22[count2] <- vector2[k]
                diff <- igraph::difference(graph, graphRewire)
                
                #Random
                Random <- rewireOnl(data=graphRandom, number=z)
                Random <- igraph::union(Random, diffR)
                comrR1 <- methodCommunity(graph=Random, method=method1,
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                comrR2 <- methodCommunity(graph=Random, method=method2,
                                         directed=directed,
                                         weights=weights,
                                         steps=steps, 
                                         spins=spins, 
                                         e.weights=e.weights, 
                                         v.weights=v.weights, 
                                         nb.trials=nb.trials)
                vectorR1[k] <- compare(comrR1, comRandom1, method="vi")
                vectorR2[k] <- compare(comrR2, comRandom2, method="vi")
                viR11[count2] <- vectorR1[k]
                viR22[count2] <- vectorR2[k]
                diffR <- igraph::difference(graphRandom, Random)
                for(k in c(2:nrep)) 
                {
                    count2 <- count2+1
                    graphRewire <- rewireOnl(data=graphRewire,
                                            number=round(0.01*z))
                    comr1 <- methodCommunity(graph=graphRewire, method=method1,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    comr2 <- methodCommunity(graph=graphRewire, method=method2,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    vector1[k] <- compare(comr1, comReal1, method="vi")
                    vector2[k] <- compare(comr2, comReal2, method="vi")
                    vi11[count2] <- vector1[k]
                    vi22[count2] <- vector2[k]
                    #Random
                    graphRewireRandom <- rewireOnl(data=Random,
                                             number=round(0.01*z))
                    comrR1 <- methodCommunity(graph=Random, method=method1,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    comrR2 <- methodCommunity(graph=Random, method=method2,
                                             directed=directed,
                                             weights=weights,
                                             steps=steps, 
                                             spins=spins, 
                                             e.weights=e.weights, 
                                             v.weights=v.weights, 
                                             nb.trials=nb.trials)
                    vectorR1[k] <- compare(comrR1, comRandom1, method="vi")
                    vectorR2[k] <- compare(comrR2, comRandom2, method="vi")
                    viR11[count2] <- vectorR1[k]
                    viR22[count2] <- vectorR2[k]
                }
                viMean11[s] <- mean(vi11)
                viMean22[s] <- mean(vi22)
                viMeanRandom11[s] <- mean(viR11)
                viMeanRandom22[s] <- mean(viR22)
            }
            graph <- igraph::intersection(graph, graphRewire)
            viMean1 <-cbind(viMean1,viMean11)
            viMean2 <-cbind(viMean2,viMean22)
            graphRandom <- igraph::intersection(graphRandom, Random)
            viMeanRandom1 <-cbind(viMeanRandom1,viMeanRandom11)
            viMeanRandom2 <-cbind(viMeanRandom2,viMeanRandom22)
            z1 <- igraph::gsize(graph)
            print(z1)
            }
    }
    colnames(viMean1) <- nRewire 
    colnames(viMean2) <- nRewire 
    colnames(viMeanRandom1) <- nRewire
    colnames(viMeanRandom2) <- nRewire
    nn <- rep(nRewire, each=nrep) 
    
    
    ratios1 <- log2((viMean1+0.001)/(viMeanRandom1+0.001))
    ratios2 <- log2((viMean2+0.001)/(viMeanRandom2+0.001))
    ratios1vs2 <- log2((viMean1+0.001)/(viMean2+0.001))
    bats1 <- as.vector(ratios1)
    bats2 <- as.vector(ratios2)
    bats1vs2 <- as.vector(ratios1vs2)
    names(bats1) <- nn
    names(bats2) <- nn
    names(bats1vs2) <- nn
    resBats1 <- cbind(ID="ratios", t(bats1))
    resBats2 <- cbind(ID="ratios", t(bats2))
    resBats1vs2 <- cbind(ID="ratios", t(bats1vs2))

    output <- list( viMean1=viMean1,
                    viMean2=viMean2,
                    viMeanRandom1=viMeanRandom1,
                    viMeanRandom2=viMeanRandom2,
                   ratios1=resBats1,
                   ratios2=resBats2,
                   ratios1vs2=resBats1vs2)
    return(output)
}


################### PLOT COMPARISON ##################
#' plotRobinCompare
#'
#' @param graph The input graph prepNet
#' @param legend The legend for the two graphs of the two methods. 
#' The default is c("real data","null model")
#' @param legend1vs2 The legend for the graph of he comparison of the two method.
#' The default is c("method1","method2")
#'
#' @return A plot
#' @export
#'
#' @examples
plotRobinCompare <- function(graph, model1, modelR1, model2, modelR2, 
                            legend=c("real data", "null model"),
                            legend1vs2=c("method1", "method2"),
                            title1="Method1",title2="Method2",
                            title1vs2="Method1 vs Method2")
{
    plot1 <- plotRobin(graph=graph, model=model1, 
                    modelR=modelR1, legend=legend)+
                    ggtitle(title1)
    plot2 <- plotRobin(graph=graph, model=model2, modelR=modelR2,
              legend=legend)+ggtitle(title2)
    plot3 <- plotRobin(graph=graph, model=model1, modelR=model2,
              legend=legend1vs2)+ggtitle(title1vs2)
    
    gridExtra::grid.arrange(plot1,plot2,plot3, nrow=2)
}

########################### TEST ROBIN ###################

################ CREATE ITPSpline ##################
#' createITPSplineResult
#' 
#' @description creates an fdatest::ITP2 class object from an iGraph
#' 
#' @param graph an iGraph object
#' @param model1 
#' @param model2 
#' @param muParam the mu parameter for ITP2bspline (default 0)
#' @param orderParam the order parameter for ITP2bspline (default 4)
#' @param nKnots the nknots parameter for ITP2bspline (default 12)
#' @param BParam the B parameter for ITP2bspline (default 10000)
#' @param isPaired the paired parameter for ITP2bspline (default TRUE)
#'
#' @return an ITP2 object
#' @export
#'
#' @examples
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

#' Title
#'
#' @param ratio 
#'
#' @return
#' @export
#'
#' @examples
robinGPTest <- function(ratio)
{
    gpregeOptions = list(indexRange=(1:2), explore=FALSE, 
                         exhaustPlotRes=30, exhaustPlotLevels=10, 
                         exhaustPlotMaxWidth=100, iters=100, 
                         labels=rep(FALSE,2), display=FALSE)
    MA=as.matrix(ratio[,c(2:dim(ratio)[2])])
    MA<-t(MA)
    vt=unique(colnames(MA))
    ntimes=length(vt)
    stdv=NULL
    varv=NULL
    for (i in c(1:ntimes)){    #ntime number of percentuage of rewire
        ind=which(colnames(MA)==vt[i])
        stdv[i]=sd(MA[1,ind])
        varv[i]=var(MA[1,ind])
    }
    #sigmaest=mean(stdv)Order of the B-spline basis expansion.
    GlobalVar=var(MA[1,])
    SigNoise=mean(varv)/GlobalVar
    if (SigNoise>1)SigNoise=1
    #
    
    #SigNoise=1-var(MA[2,])
    sigmaest=1-SigNoise
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
    
    dvet=data.matrix(as.numeric(colnames(ratio)[-1]))
    dd=t(data.matrix(as.numeric((ratio)[-1])))
    rownames(dd)='VI'
    colnames(dd)=dvet
    datadum=rbind(dd,dd)
    gpregeOutput <-gprege::gprege(data=datadum, inputs=dvet,
                                  gpregeOptions=gpregeOptions)
    
    bf=gpregeOutput$rankingScores[1]
    return(Bayes_Factor=bf)
 }


################ FDA TEST ##############################

#' robinFDATest
#'
#' @param graph The input graph created with prepGraph
#' @param model1 
#' @param model2 
#' @param legend 
#'
#' @return
#' @export
#'
#' @examples
robinFDATest <- function(graph,
                        model1=model1,
                        model2=model2, 
                        legend=legend)
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
        ggplot2::ylab("Variation of Information (VI)")+
        ggplot2::ggtitle("Robin plot")
    plot1
    print(plot1)
   
    perc<-rep((seq(0,60,5)/100))
    ITPresult <- createITPSplineResult(graph, model1, model2)
    plot2<-plot(ITPresult, main='VI', xrange=c(0,0.6), xlab='perturbation', 
                ylab="VI")
    lines(perc, rep(0.05, 13), type="l", col="red")
    
    print(plot2)
}  

########### AREA UNDER THE CURVE    ##############
#' Title
#'
#' @param graph The input graph created with prepGraph
#' @param model1 
#' @param model2 
#'
#' @return
#' @export
#'
#' @examples
robinAUCTest <- function(graph,
                         model1=model1,
                         model2=model2)
{
    N <- igraph::vcount(graph)
    mvimeanmodel1 <- cbind(as.vector((apply(model1, 2, mean))/log2(N)))
    mvimeanmodel2 <- cbind(as.vector((apply(model2, 2, mean))/log2(N)))
    area1<-DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel1, 
                          method ="spline")
    area2<-DescTools::AUC(x=(seq(0,60,5)/100), y=mvimeanmodel2, 
                          method ="spline")
    output <- list(area1=area1,area2=area2)
return(output)
}





    