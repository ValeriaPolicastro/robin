######PREPARAZIONE NETWORK########## 
prepNet <- function(net){
    edge<-read.table(net, quote="\"")
    edge<-as.matrix(edge)
    ee<-edge
    vet1<-as.vector(t(edge))
    un<-unique(sort(vet1))
    if(max(un)!=length(un)){ ##se non vi sono tutti i nodi il massimo del vettore non è uguale alla lunghezza
        id<-seq(1,length(un))  ##il nome dei nodi è vet che è uguale all'id, altrimenti è vet1
        vet<-vet1
        for(i in c(1:length(un))){
            ind<-which(vet1==un[i])
            vet[ind]<-id[i]
            }
        edge<-matrix(vet,ncol=2,byrow=TRUE)
        graph<-graph(vet, directed=FALSE )
        }else{
            graph<-graph(vet1, directed=FALSE )
            }
    graph<-simplify(graph) #grafici che non contengono loop ed archi multipli
    return(graph)
}
#tenere gli id?
 
 ######MODELLO RANDOM#########
random<-function(graph){
   z<- ecount(graph) #numero di archi
   graphRandom<-rewire(graph,with=keeping_degseq(loops = FALSE, niter = z))
   #Riordinare casualmente i vertici preservando i gradi del grafico
   return(graphRandom)
 }

#####TUTTI I METODI DI IGRAPH####    
methodCommunity<- function(graph,type,weights = NULL,steps = 4,spins = 25,e.weights = NULL,v.weights = NULL, nb.trials = 10,directed = TRUE){
    switch(type,
           louvain=cluster_louvain(graph=graph,weights = NULL),
           walktrap=cluster_walktrap(graph=graph,
                                     weights = E(graph)$weight,steps = 4),
           spinglass=cluster_spinglass(graph=graph, weights = NULL,spins = 25),
           leadingEigen=cluster_leading_eigen(graph=graph, steps = -1,
                                              weights = NULL),
           edgeBetweenness=cluster_edge_betweenness(graph=graph,
                                                    weights = E(graph)$weight, directed = TRUE),
           fastGreedy=cluster_fast_greedy(graph=graph,weights = E(graph)$weight),
           #usato questo loro fastgreedy.community(graph)
           labelProp=cluster_label_prop(graph=graph, weights = NULL),
           infomap=cluster_infomap(graph=graph, e.weights = NULL,
                                   v.weights = NULL, nb.trials = 10)
    )
}
com<-membership(methodCommunity(graph=graph,type="louvain"))
# com2<-membership(methodCommunity(graph,type="walktrap"))
# com3<-membership(methodCommunity(graph,type="spinglass"))
# com4<-membership(methodCommunity(graph,type="leadingEigen",steps =-1))
# com5<-membership(methodCommunity(graph,type="edgeBetweenness"))
# com6<-membership(methodCommunity(graph,type="fastGreedy"))
# com7<-membership(methodCommunity(graph,type="labelProp"))
# com8<-membership(methodCommunity(graph,type="infomap"))
# 
# com1
# com2
# com3
# com4
# com5
# com6
# com7
# com8
##non sono tutti i parametri dati e la nostra è directed FALSE se gli metto nello switch i parametri di defoult se li prende? 




########OUTPUT#######
iter <- function(base, graph, graphRandom, methodCommunity) {
    fileoutbats <- paste(base,"_BATS.txt",sep="")
    filepdf <- paste(base,"_VI.pdf",sep="")
    fileoutvirand <- paste(base,"_VI_random.txt",sep="")
    fileoutvicase <- paste(base,"_VI_case.txt",sep="")
    #file utilizzati da Italia per FAD
    fileoutvirandbio <- paste(base,"_VI_random_bio.txt",sep="")
    fileoutvicasebio <- paste(base,"_VI_case_bio.txt",sep="")
    nrep <- 10
    N <- vcount(graph)  ##numero di vertici


####RETE VERA#
if(method=="louvain"){  #cercare le community con louvain o con fastgreedy
  ebc <- cluster_louvain(graph)
  }else{
      ebc <- fastgreedy.community(graph)
      }
    com<-membership(ebc)#gives the division of the vertices, into communities

####RETE RANDOM#
if(method=="louvain"){ 
   ebc <- cluster_louvain(graphRandom)
}else{
  ebc <- fastgreedy.community(graphRandom)
  }
    comn <- membership(ebc)
    de <- ecount(graph)

vet1 <- seq(5,100,5) #dal 5 a 100 con passo 5 
vet <- round(vet1*de/100,0)#arrotonda a 0 cifre decimali

#MATRICI OUTPUT
vi_bhl <- matrix(0,nrep^2,length(vet)+1)
vi_random <- matrix(0,nrep^2,length(vet)+1)
vi_mean_random <- matrix(0,nrep,length(vet)+1)
vi_mean_bhl <- matrix(0,nrep,length(vet)+1)

count <- 1
for(z in vet){
     count2 <- 0
     vi <- NULL
     count <- count+1
     for(s in c(1:nrep)){#ciclo da 1 a 10 
         count2 <- count2+1
         vet_random <- NULL
     vet_bhl <- NULL
     k <- 1
     ###perturba il grafo reale e ricalcola le community
     g3 <- rewire(graph,with=keeping_degseq(loops = FALSE, niter = z))
     
     if(method=="louvain"){
       ebc <- igraph::cluster_louvain(g3)
     }else{
       ebc <- igraph::fastgreedy.community(g3)
     }
     comr <- membership(ebc)
     #compare the community of real and rewired through VI
     vet_bhl[k] <- compare(com,comr,method="vi")
     vi_bhl[count2,count] <- vet_bhl[k]
     ###perturba il grafo random e ricalcola le community
     g3_random <- rewire(graphRandom,with=keeping_degseq(loops = FALSE, niter = z))
     if(method=="louvain"){
        ebc <- cluster_louvain(g3_random)
     }else{
       ebc <- fastgreedy.community(g3_random)
     }
     comr <- membership(ebc)
     #compare the community made with the random model and the rewired
     vet_random[k] <- compare(comn,comr,method="vi")
     vi_random[count2,count] <- vet_random[k]
     for(k in c(2:nrep)){ #lo itera 10 volte non 100
         count2 <- count2+1
         g3_iter <- rewire(g3,with=keeping_degseq(loops = FALSE, niter = round(0.01*z)))
     if(method=="louvain"){
        ebc <- cluster_louvain(g3_iter)
     }else{
        ebc <- fastgreedy.community(g3_iter)
     }
     comr <- membership(ebc)
     #compare distance real and rewired through VI
     vet_bhl[k] <- compare(com,comr,method="vi")
     vi_bhl[count2,count] <- vet_bhl[k]
     g3_rand_it <- rewire(g3_random,with=keeping_degseq(loops = FALSE, niter = round(0.01*z)))
     if(method=="louvain"){
         ebc <- cluster_louvain(g3_rand_it)
     }else{
        ebc <- fastgreedy.community(g3_rand_it)
     }   
     comr <- membership(ebc)
     #compare random and rewired random through VI
     vet_random[k] <- compare(comn,comr,method="vi")
     vi_random[count2,count] <- vet_random[k]
     }
    vi_mean_random[s,count] <- mean(vet_random)
    vi_mean_bhl[s,count] <- mean(vet_bhl)
     }
     }
nn1 <- c(0,vet1)

colnames(vi_random) <- nn1
colnames(vi_bhl) <- nn1

write.table(vi_random,fileoutvirand,sep="\t",row.names=FALSE,quote=FALSE)
write.table(vi_bhl,fileoutvicase,sep="\t",row.names=FALSE,quote=FALSE)

colnames(vi_mean_random) <- nn1
colnames(vi_mean_bhl) <- nn1

write.table(vi_mean_random,fileoutvirandbio,sep="\t",row.names=FALSE,quote=FALSE)
write.table(vi_mean_bhl,fileoutvicasebio,sep="\t",row.names=FALSE,quote=FALSE)

#rapporto tra la media delle distanze VI tra il modello reale e quello perturbato e la media delle distanze tra il random e la sua perturbazione
ratios <- log2((vi_mean_bhl+0.001)/(vi_mean_random+0.001))
bats <- as.vector(ratios)
nn <- rep(c(0,vet1),each=nrep)
names(bats) <- nn

resBats <- cbind(ID="ratios",t(bats))#la trasposta del rapporto
write.table(resBats,fileoutbats,sep="\t",row.names=FALSE,quote=FALSE)
mvi_bhl <- apply(vi_mean_bhl,2,mean)
mvi_random <- apply(vi_mean_random,2,mean)
vv <- c(0,vet1/100)

pdf(filepdf)

####PLOT
plot(vv,mvi_random/log2(N),col="red",type="o",axes=FALSE, ann=FALSE,ylim=c(0,1))
axis(1,at=c(0,0.2,0.4,0.6,0.8,1),lab=c(0,0.2,0.4,0.6,0.8,1))
axis(2,las=1,at=c(0,0.2,0.4,0.6,0.8,1),lab=c(0,0.2,0.4,0.6,0.8,1))
box()
#abline(h=c(0.1,0.2))
lines(vv,mvi_bhl/log2(N),type="o",pch=22,lty=2,col="blue")
legend("topleft",pch=c("o","-"),legend=c("null model","real data"),col=c("red","blue"))
title(xlab="percentage of perturbation")
title(ylab="variation of information (VI)")
dev.off()
}
