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
        id <- seq(1,length(un))  ##il nome dei nodi è vet che è uguale all'id, altrimenti è vet1
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
 
######MODELLO RANDOM#########
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
    z <- ecount(graph) #number of edges
    graphRandom <- rewire(graph, with=keeping_degseq(loops=FALSE, niter=z))
    #Randomly rewire the edges while preserving the original graph's degree distribution
    #for z all the edges
    return(graphRandom)
}

#####COMMUNITY METHOD####    
#' methodCommunity
#' @description This function gives the membership vector of the community structure. The community structure was found by functions implemented in igraph.
#' @param graph the input graph
#' @param method the clustering method, one of "walktrap", "edgeBetweenness", 
#' "fastGreedy", "louvain", "spinglass", "leadingEigen","labelProp","infomap"
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
                            directed=TRUE){

    if(is.null(weights) &
        (sum(method %in% c("walktrap", "edgeBetweenness", "fastGreedy") == 1 )))
    {
        weights <- E(graph)$weight
    }

    if((steps == 4) & (method == "leadingEigen"))
     {
          steps <--1
      }
    communities <- switch(method,
           louvain=cluster_louvain(graph=graph,weights=weights),
           walktrap=cluster_walktrap(graph=graph,weights=weights,steps=steps),
           spinglass=cluster_spinglass(graph=graph,weights=weights,spins=spins),
           leadingEigen=cluster_leading_eigen(graph=graph,steps=steps,
                                              weights=weights),
           edgeBetweenness=cluster_edge_betweenness(graph=graph,
                                                    weights=weights,
                                                    directed=directed),
           fastGreedy=cluster_fast_greedy(graph=graph,weights=weights),
           labelProp=cluster_label_prop(graph=graph,weights=weights),
           infomap=cluster_infomap(graph=graph,e.weights=e.weights,
                                   v.weights=v.weights,nb.trials=nb.trials)
    )
    return(membership(communities))
}

##la nostra è directed FALSE se gli metto nello switch i parametri di defoult se li prende? 


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
iter<- function(base, graph, graphRandom,method,type) {
    nrep <- 10
    N <- vcount(graph) 
    com<-methodCommunity(graph,method)#real network
    comn<-methodCommunity(graphRandom,method)#random network
    de<-ecount(graph)
    vet1<-seq(5,100,5) #dal 5 a 100 con passo 5 
    vet<-round(vet1*de/100,0)#arrotonda a 0 cifre decimali
 #OUTPUT MATRIX
    viBhl<-matrix(0,nrep^2,length(vet)+1)
    viRandom<-matrix(0,nrep^2,length(vet)+1)
    viMeanRandom<-matrix(0,nrep,length(vet)+1)
    viMeanBhl<-matrix(0,nrep,length(vet)+1)
    
    count<-1
 #INDEPENDENT    
    if(type=="independent"){ 
        for(z in vet){
            count2<-0
            vi<-NULL
            count<-count+1
            for(s in c(1:nrep)){
                count2<-count2+1
                vetRandom<-NULL
                vetBhl<-NULL
                k<-1
        #perturba il grafo reale e ricalcola le community
                graphRewire<-rewire(graph,
                                    with=keeping_degseq(loops = FALSE,
                                                        niter = z))
                comr <-methodCommunity(graphRewire,method)
                #compare the community of real and rewired through VI
                vetBhl[k]<-compare(com,comr,method="vi")
            viBhl[count2,count]<-vetBhl[k]
        #perturba il grafo random e ricalcola le community
                graphRewireRandom<-rewire(graphRandom,
                            with=keeping_degseq(loops = FALSE, niter = z))
                comr<- methodCommunity(graphRewireRandom,method)
                #compare the community made with the random 
                #model and the rewired
                vetRandom[k]<-compare(comn,comr,method="vi")
            viRandom[count2,count]<-vetRandom[k]
                for(k in c(2:nrep)){
                    count2<-count2+1
                    graphRewireIter<-rewire(graphRewire,
                    with=keeping_degseq(loops = FALSE, niter = round(0.01*z)))
                    comr<-methodCommunity(graphRewireIter,method)
                    #compare distance real and rewired through VI
                    vetBhl[k]<-compare(com,comr,method="vi")
                viBhl[count2,count]<-vetBhl[k]
                    graphRewireRandIt<-rewire(graphRewireRandom,
                    with=keeping_degseq(loops = FALSE, niter = round(0.01*z)))
                    comr <- methodCommunity(graphRewireRandIt,method)
                    #compare random and rewired random through VI
                    vetRandom[k]<-compare(comn,comr,method="vi")
                viRandom[count2,count]<-vetRandom[k]
                }
            viMeanRandom[s,count]<-mean(vetRandom)
            viMeanBhl[s,count]<-mean(vetBhl)
            }
        }
#DEPENDENT 
    }else{
        (type=="dependent")
    vet<-round((5*ecount(graph))/100,0)
    diff<-NULL
    diffR<-NULL
    diffRI<-NULL
        for(z in vet){
                        count2<-0
                        vi<-NULL
                        count<-count+1
                   #for(s in c(1:nrep)){
                            count2<-count2+1
                            vetRandom<-NULL
                            vetBhl<-NULL
                           #k<-1
                        ###REAL
                            #perturba il grafo reale e ricalcola le community
                            graphRewire<-rewire(graph,
                                                with=keeping_degseq(loops = FALSE,
                                                                    niter = z))
                            graphRewire<-union(graphRewire,diff)
                            comr<-methodCommunity(graphRewire,method)
                            vetBhl[k]<-compare(com,comr,method="vi")
                            viBhl[count2,count]<-vetBhl[k]
                            diff<-difference(graph,graphRewire)
                            graph<-intersection(graph,graphRewire)
                            vet<-round((5*ecount(graph))/100,0)
                            com<-comr
                            
                            
                        ###RANDOM    
                            # graphRewireRandom<-rewire(graphRandom,
                            #                           with=keeping_degseq(loops = FALSE, niter = z))
                            # graphRewireRandom<-union(graphRewire,diffR)
                            # comr<- methodCommunity(graphRewireRandom,method)
                            # vetRandom[k]<-compare(comn,comr,method="vi")
                            # viRandom[count2,count]<-vetRandom
                            # diffR<-difference(graphRandom,graphRewireRandom)
                            # graphRandom<-intersection(graphRandom,graphRewireRandom)
                            # vetR<-round((5*ecount(graphRandom))/100,0)
                            # #deve avere un for da solo questo vet è diverso dall' altro 
                            # comn<-comr
                            
     

                            # for(k in c(2:nrep)){
                            #     count2<-count2+1
                            #     graphRewireIter<-rewire(graphRewire,
                            #        with=keeping_degseq(loops = FALSE, niter = round(0.01*z)))
                            #     graphRewire<-union(graphRewireInter,diffRI)
                            #     comr<-methodCommunity(graphRewireIter,method)
                            #     #compare distance real and rewired through VI
                            #     vetBhl[k]<-compare(com,comr,method="vi")
                            #     viBhl[count2,count]<-vetBhl[k]
                            #     diffRI<-difference(graphRewire,graphRewireIter)
                            #     graph<-intersection(graphRewire,graphRewireIter)
                            #     vetIR<-round((5*ecount(graph))/100,0)
                            #     com<-comr
                            #     
                            #######RANDOM    
                            #     graphRewireRandIt<-rewire(graphRewireRandom,
                            #      with=keeping_degseq(loops = FALSE, niter = round(0.01*z)))
                            #     comr <- methodCommunity(graphRewireRandIt,method)
                            # vetRandom[k]<-compare(comn,comr,method="vi")
                            #     viRandom[count2,count]<-vetRandom[k]
                            }
                           # viMeanRandom[s,count]<-mean(vetRandom)
                            viMeanBhl[s,count]<-mean(vetBhl)
                   }
        }
    
    
         
    
    
    
    fileoutbats<-paste(base,"_BATS.txt",sep="")
    filepdf<-paste(base,"_VI.pdf",sep="")
    fileoutvirand<-paste(base,"_VI_random.txt",sep="")
    fileoutvicase<-paste(base,"_VI_case.txt",sep="")
    #file utilizzati da Italia per FAD
    fileoutvirandbio<-paste(base,"_VI_random_bio.txt",sep="")
    fileoutvicasebio<-paste(base,"_VI_case_bio.txt",sep="")
    
    nn1<-c(0,vet1)
    colnames(viRandom) <- nn1
    colnames(viBhl) <- nn1
    write.table(viRandom,fileoutvirand,sep="\t",row.names=FALSE,quote=FALSE)
    write.table(viBhl,fileoutvicase,sep="\t",row.names=FALSE,quote=FALSE)
    colnames(viMeanRandom) <- nn1
    colnames(viMeanBhl) <- nn1
    write.table(viMeanRandom,fileoutvirandbio,sep="\t",row.names=FALSE,quote=FALSE)
    write.table(viMeanBhl,fileoutvicasebio,sep="\t",row.names=FALSE,quote=FALSE)
    
    #rapporto tra la media delle distanze VI tra il modello reale e quello
    ##perturbato e la media delle distanze tra il random e la sua perturbazione
    ratios<-log2((viMeanBhl+0.001)/(viMeanRandom+0.001))
    bats<-as.vector(ratios)
    nn<-rep(c(0,vet1),each=nrep)
    names(bats)<-nn
    resBats<-cbind(ID="ratios",t(bats))#la trasposta del rapporto
    write.table(resBats,fileoutbats,sep="\t",row.names=FALSE,quote=FALSE)
    mviBhl <- apply(viMeanBhl,2,mean)
    mviRandom <- apply(viMeanRandom,2,mean)
   
####PLOT#######
    vv<-c(0,vet1/100)
    pdf(filepdf)
    plot(vv,mviRandom/log2(N),col="red",type="o",axes=FALSE, ann=FALSE,ylim=c(0,1))
    axis(1,at=c(0,0.2,0.4,0.6,0.8,1),lab=c(0,0.2,0.4,0.6,0.8,1))
    axis(2,las=1,at=c(0,0.2,0.4,0.6,0.8,1),lab=c(0,0.2,0.4,0.6,0.8,1))
    box()
    #abline(h=c(0.1,0.2))
    lines(vv,mviBhl/log2(N),type="o",pch=22,lty=2,col="blue")
    legend("topleft",pch=c("o","-"),legend=c("null model","real data"),col=c("red","blue"))
    title(xlab="percentage of perturbation")
    title(ylab="variation of information (VI)")
    dev.off()
}
###COMPARISON DIFFERENT METHODS####
comparison<-function(method1,method2){
    fist<-iter()
    second<-
        }
