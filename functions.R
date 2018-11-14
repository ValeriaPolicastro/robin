######PREPARAZIONE NETWORK########## 
#' Title
#'
#' @param net network
#'
#' @return simple graph
#' @export
#'
#' @examples
prepNet <- function(net) 
{
    edge <- read.table(net, quote="\"")
    edge <- as.matrix(edge)
    ee <- edge
    vet1 <- as.vector(t(edge))
    un <- unique(sort(vet1))
    if(max(un) != length(un)) { ##se non vi sono tutti i nodi il massimo del vettore non è uguale alla lunghezza
        id <- seq(1, length(un))  ##il nome dei nodi è vet che è uguale all'id, altrimenti è vet1
        vet <- vet1
        for(i in c(1:length(un))) {
            ind <- which(vet1 == un[i])
            vet[ind] <- id[i]
        }
        edge <- matrix(vet, ncol=2, byrow=TRUE)
        graph <- graph(vet, directed=FALSE)
    }else{
        graph <- graph(vet1, directed=FALSE)
    }
    graph <- simplify(graph) #grafici che non contengono loop ed archi multipli
    return(graph)
}
#tenere gli id?
 
######GRAPH RANDOM#########
#' graphRandom
#'
#' @param graph graph data
#'
#' @return graph random
#' @export
#'
#' @examples
random <- function(graph)
{
    z <- ecount(graph) ## number of edges
    graphRandom <- rewire(graph, with=keeping_degseq(loops=FALSE, niter=z))
    ## Randomly rewire the edges while preserving the original graph's degree distribution
    ## for z all the edges
    return(graphRandom)
}

#####COMMUNITY METHOD####    
#' methodCommunity
#' @description This function gives the membership vector of the community structure. The community structure was found by functions implemented in igraph.
#' @param graph the input graph
#' @param method the clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen", "labelProp", "infomap"
#' @param weights this argument is not settable for "infomap" method
#' @param steps this argument is settable only for "leadingEigen"and"walktrap" method
#' @param spins this argument is settable only for "infomap" method
#' @param e.weights this argument is settable only for "infomap" method
#' @param v.weights this argument is settable only for "infomap" method
#' @param nb.trials this argument is settable only for "infomap" method
#' @param directed this argument is settable only for "edgeBetweenness" method
#'
#' @return membership vector of the community structure
#' @export
#'
#' @examples
methodCommunity <- function(graph, 
                            method, 
                            weights=NULL, 
                            steps=4, 
                            spins=25, 
                            e.weights=NULL, 
                            v.weights=NULL, 
                            nb.trials=10, 
                            directed=TRUE)
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
           louvain=cluster_louvain(graph=graph, 
                                    weights=weights), 
           walktrap=cluster_walktrap(graph=graph, 
                                    weights=weights, 
                                    steps=steps), 
           spinglass=cluster_spinglass(graph=graph, 
                                        weights=weights, 
                                        spins=spins), 
           leadingEigen=cluster_leading_eigen(graph=graph, 
                                            steps=steps, 
                                            weights=weights), 
           edgeBetweenness=cluster_edge_betweenness(graph=graph, 
                                                    weights=weights, 
                                                    directed=directed), 
           fastGreedy=cluster_fast_greedy(graph=graph, weights=weights), 
           labelProp=cluster_label_prop(graph=graph, weights=weights), 
           infomap=cluster_infomap(graph=graph, e.weights=e.weights, 
                                v.weights=v.weights, nb.trials=nb.trials)
    )
    return(membership(communities))
}

##la nostra è directed FALSE se gli metto nello switch i parametri di defoult se li prende? 

#########REWIRE COMPLETE
#' Title
#'
#' @param data 
#' @param number 
#' @param community 
#'
#' @return
#' @export
#'
#' @examples
rewireCompl <- function(data, number, community, method)
{
    graphRewire <- rewire(data, with=keeping_degseq(loops=FALSE, niter=number))
    comR  <- methodCommunity(graph=graphRewire, method=method)
    VI <- compare(community, comR, method="vi")
    output <- list(VI=VI, graphRewire=graphRewire)
    return(output)#non serve per tutti
}
#perturba il grafo e ricalcola le community
#compare the community through VI

#########REWIRE ONLY
#' Title
#'
#' @param data 
#' @param number 
#'
#' @return
#' @export
#'
#' @examples
rewireOnl <- function(data, number)
{
    graphRewire <- rewire(data, with=keeping_degseq(loops = FALSE, niter = number))
    return(graphRewire)
}

############PLOT##############
#' Title
#'
#' @param graph 
#'
#' @return
#' @export
#'
#' @examples
plotRobin <-  function(graph, base="indip")
{   
    filepdf <- paste(base, "_VI.pdf", sep="")
    viMeanBhl <- iter(base="independent", graph=graph, graphRandom=graphRandom, 
                    method="fastGreedy", type="independent")[[3]]
    viMeanRandom <- iter(base="independent", graph=graph, graphRandom=graphRandom, 
                       method="fastGreedy", type="independent")[[4]]
    N <- vcount(graph)
    mviBhl <- apply(viMeanBhl, 2, mean)
    mviRandom <- apply(viMeanRandom, 2, mean)
    vet1 <- seq(5, 100, 5)
    vv <- c(0, vet1/100)
    plotVI <- pdf(filepdf)
    plot(vv, mviRandom/log2(N), col="red", type="o", axes=FALSE, ann=FALSE, ylim=c(0, 1))
    axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), lab=c(0, 0.2, 0.4, 0.6, 0.8, 1))
    axis(2, las=1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), lab=c(0, 0.2, 0.4, 0.6, 0.8, 1))
    box()
    #abline(h=c(0.1, 0.2))
    lines(vv, mviBhl/log2(N), type="o", pch=22, lty=2, col="blue")
    legend("topleft", pch=c("o", "-"), legend=c("null model", "real data"), col=c("red", "blue"))
    title(xlab="percentage of perturbation")
    title(ylab="variation of information (VI)")
    dev.off() 
    return(plotVI)
}

########ITER#######
#' Title
#'
#' @param base name of the output
#' @param graph graph data
#' @param graphRandom graph random data
#' @param method community method
#' @param type dependent or independent  
#'
#' @return
#' @export
#'
#' @examples
iter <- function(base, graph, graphRandom, method, 
                type=c("dependent", "independent")) 
{
    type <- match.arg(type)
    nrep <- 10
    comReal <- methodCommunity(graph=graph, method=method) # real network
    comRandom <- methodCommunity(graph=graphRandom, method=method) # random network
    de <- ecount(graph)
    VI <- NULL
    vetBhl <- NULL
    vetRandom <- NULL
    graphRewireRandom <- NULL
    graphRewire <- NULL
    count <- 1
    vet1 <- seq(5, 100, 5)  #dal 5 a 100 con passo 5 
 #INDEPENDENT    
    if(type == "independent") 
    {
        vet <- round(vet1*de/100, 0)#arrotonda a 0 cifre decimali
        #OUTPUT MATRIX
        viBhl <- matrix(0, nrep^2, length(vet)+1)
        viRandom <- matrix(0, nrep^2, length(vet)+1)
        viMeanRandom <- matrix(0, nrep, length(vet)+1)
        viMeanBhl <- matrix(0, nrep, length(vet)+1)
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
                                    method=method)
                vetBhl[k] <- Real$VI
                graphRewire <- Real$graphRewire
                viBhl[count2, count] <- vetBhl[k]
                #RANDOM
                Random <- rewireCompl(data=graphRandom, number=z, 
                                  community=comRandom, method=method)
                vetRandom[k] <- Random$VI
                viRandom[count2, count] <- vetRandom[k]
                graphRewireRandom <- Random$graphRewire
                for(k in c(2:nrep))
                {
                    count2 <- count2+1
                    Real <- rewireCompl(data=graphRewire, 
                                    number=round(0.01*z), 
                                    community=comReal, 
                                    method=method)
                    vetBhl[k] <- Real$VI
                    viBhl[count2, count] <- vetBhl[k]
                    Random <- rewireCompl(data=graphRewireRandom, 
                                     number=round(0.01*z), community=comRandom, 
                                     method=method)
                    vetRandom[k] <- Random$VI
                    viRandom[count2, count] <- vetRandom[k]
                }
                viMeanRandom[s, count] <- mean(vetRandom)
                viMeanBhl[s, count] <- mean(vetBhl)
            }
        }
#DEPENDENT 
    }else{
        
        vet1 <- round((5*ecount(graph))/100, 0)
        vet <- rep(vet1, 19) #19 volte più quello che manca all' ultimo giro
        #ultimo <- de-(sum(vet))
        viBhl <- matrix(0, nrep^2, length(vet1)+1)
        viRandom <- matrix(0, nrep^2, length(vet1)+1)
        viMeanRandom <- matrix(0, nrep, length(vet1)+1)
        viMeanBhl <- matrix(0, nrep, length(vet1)+1)
        diff <- NULL
        diffR <- NULL
        for(z in vet)
        {
            count2 <- 0
            count <- count+1
            for(s in c(1:nrep)){
                count2 <- count2+1
                k <- 1
                ###REAL
                graphRewire <- rewireOnl(data=graph, number=z)
                graphRewire <- union(graphRewire, diff)
                comr <- methodCommunity(graph=graphRewire, method=method)
                vetBhl[k] <- compare(comReal, comr, method="vi")
                viBhl[count2, count] <- vetBhl[k]
                diff <- difference(graph, graphRewire)
                graph <- intersection(graph, graphRewire)
                z1 <- round((5*ecount(graph))/100, 0)
                ###RANDOM  
                graphRewireRandom <- rewireOnl(data=graphRandom, number=z)
                graphRewireRandom <- union(graphRewireRandom, diffR)
                comr <- methodCommunity(graph=graphRewireRandom, method=method)
                vetBhl[k] <- compare(comRandom, comr, method="vi")
                viBhl[count2, count] <- vetBhl[k]
                diffR <- difference(graphRandom, graphRewireRandom)
                graphRandom <- intersection(graphRandom, graphRewireRandom)
                z2 <- round((5*ecount(graphRandom))/100, 0)
                for(k in c(2:nrep)) 
                {
                    count2 <- count2+1
                    #REAL
                    Real <- rewireCompl(data=graphRewire, number=round(0.01*z),
                                        method=method,
                                   community=comReal)
                    vetBhl[k] <- Real$VI
                    viBhl[count2, count] <- vetBhl[k]
                    ## RANDOM
                    Random <- rewireCompl(data=graphRewireRandom,
                                          method=method,
                                      number=round(0.01*z), community=comRandom)
                    vetRandom[k] <- Random$VI
                    viRandom[count2, count] <- vetRandom[k]
                }
                viMeanRandom[s, count] <- mean(vetRandom)
                viMeanBhl[s, count] <- mean(vetBhl)
            }
            }
        }
        nn1 <- c(0, vet1)
        colnames(viRandom) <- nn1
        colnames(viBhl) <- nn1
        colnames(viMeanRandom) <- nn1
        colnames(viMeanBhl) <- nn1
       
       # fileoutbats <- paste(base, "_BATS.txt", sep="")
       #  filepdf <- paste(base, "_VI.pdf", sep="")
       #  fileoutvirand <- paste(base, "_VI_random.txt", sep="")
       #  fileoutvicase <- paste(base, "_VI_case.txt", sep="")
       #  #file utilizzati da Italia per FAD
       #  fileoutvirandbio <- paste(base, "_VI_random_bio.txt", sep="")
       #  fileoutvicasebio <- paste(base, "_VI_case_bio.txt", sep="")
      
      
       #  write.table(viRandom, fileoutvirand, sep="\t", row.names=FALSE, quote=FALSE)
       #  write.table(viBhl, fileoutvicase, sep="\t", row.names=FALSE, quote=FALSE)
       #  write.table(viMeanRandom, fileoutvirandbio, sep="\t", row.names=FALSE, quote=FALSE)
       #  write.table(viMeanBhl, fileoutvicasebio, sep="\t", row.names=FALSE, quote=FALSE)
       #pdf(filepdf)
    
        ratios <- log2((viMeanBhl+0.001)/(viMeanRandom+0.001))
        #rapporto tra la media delle distanze VI tra il modello reale e quello
        ##perturbato e la media delle distanze tra il random e la sua perturbazione
        bats <- as.vector(ratios)
        nn <- rep(c(0, vet1), each=nrep)
        names(bats) <- nn
        resBats <- cbind(ID="ratios", t(bats))#la trasposta del rapporto
        #write.table(resBats, fileoutbats, sep="\t", row.names=FALSE, quote=FALSE)
        output <- list(viBhl=viBhl, viRandom=viRandom, viMeanBhl=viMeanBhl, viMeanRandom=viMeanRandom, resBats=resBats) 
        return(output)

}

###COMPARISON DIFFERENT METHODS####
####con la stessa perturbazione bisogna testare due metodi
####dipendente o indipendente??
#' Title
#'
#' @param data 
#' @param method1 
#' @param method2 
#'
#' @return
#' @export
#'
#' @examples
 comparison <- function(graph, method1, method2){
    nrep <- 10
    de <- ecount(graph)
    vet1 <- round((5*de)/100, 0)
    vet <- rep(vet1, 19)
    viBhl <- matrix(0, nrep^2, length(vet)+1)
    viMeanBhl <- matrix(0, nrep, length(vet)+1)
    VI <- NULL
    vetBhl <- NULL
    graphRewire <- NULL
    count <- 1
    diff <- NULL
    for(s in c(1:nrep)){
        count2 <- 0
        vi <- NULL
        count <- count+1
        for(z in vet){
            count2 <- count2+1
            k <- 1
            graphRewire <- rewireOnl(data=graph, number=z)
            graphRewire <- union(graphRewire, diff)
            comr1 <- methodCommunity(graph=graphRewire, method=method1)
            comr2 <- methodCommunity(graph=graphRewire, method=method2)
            vetBhl[k] <- compare(comr1, comr2, method="vi")
            viBhl[count2, count] <- vetBhl[k]
            diff <- difference(graph, graphRewire)
            graph <- intersection(graph, graphRewire)
            for(k in c(2:nrep)){
                count2 <- count2+1
                graphRewire <- rewireOnl(data=graphRewire, number=z)
                comr1 <- methodCommunity(graph=graphRewire, method=method1)
                comr2 <- methodCommunity(graph=graphRewire, method=method2)
                vetBhl[k] <- compare(comr1, comr2, method="vi")
                viBhl[count2, count] <- vetBhl[k]
            }
            viMeanBhl[s, count] <- mean(vetBhl)
        }
    }
    return(viMeanBhl)
 }
