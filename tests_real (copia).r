library("refund")
library("fdatest")
library("aqfig")


  random_bio <- read.delim(file="/media/vpoli/MYFILES/CNR Tigem/robin vecchio e alto/output_VI_random_bio.txt", header=FALSE)
  case_bio <- read.delim(file="/media/vpoli/MYFILES/CNR Tigem/robin vecchio e alto/output_VI_case_bio.txt", header=FALSE)
  
  file_out_plot=paste(dir_input,"fda_test_",i,"_",method,".eps",sep="")
  file_out_plot1=paste(dir_input,"fda_test_",i,"_",method,"_VI.eps",sep="")
  file_out_plot2=paste(dir_input,"fda_test_",i,"_",method,"_VI_functional.eps",sep="")
  file_out_plot3=paste(dir_input,"fda_test_",i,"_",method,"_pvalues.eps",sep="")
  
  file_out_test=paste(dir_input,"FAD_test_",i,"_",method,".txt",sep="")

  ntimes=20
  nn<-log2(igraph::gorder(graph))
  case_bio=as.matrix(case_bio[-1,])/nn
  random_bio=as.matrix(random_bio[-1,])/nn
 
  
  #-----------------FDATEST (Pini-Vantini paper)----------
  
  ITP.result=ITP2bspline(case_bio,random_bio,mu=0,order=4,
                         nknots=12,B=10000,paired=TRUE)
  #  pvalues_ad=ITP.result$corrected.pval
  pvalues=ITP.result$pval
  pvalues_ad=p.adjust(pvalues, method = "bonferroni", n = length(pvalues))
  
  
  # par(mfrow=c(2,2))
  # matplot(c(0:ntimes),t(random_bio),type="l", col="red",xlab='perturbation',ylab="VI",main="VI : Sampled Data")
  # matlines(c(0:ntimes),t(case_bio),type="l", col="black",xlab='perturbation',ylab="VI",main="VI : Sampled Data")
  # plot(ITP.result,main='VI',xrange=c(0,ntimes),xlab='perturbation',ylab="VI")
  
  abscissa.range=c(0,1)
  min.ascissa=abscissa.range[1]-(abscissa.range[2]-abscissa.range[1])/2
  max.ascissa=abscissa.range[2]+(abscissa.range[2]-abscissa.range[1])/2
  p=dim(ITP.result$heatmap.matrix)[1]
  ordinata.grafico <- seq(abscissa.range[1],abscissa.range[2],length.out=p) - abscissa.range[1]
  nlevel=20
  colori=rainbow(nlevel,start=0.15,end=0.67)
  colori=colori[length(colori):1]
  
  #pdf(file_out_plot)
  setEPS()
  postscript(file_out_plot)
  
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),widths=c(1,1,1))
  
  xx_ntimes=c(0:ntimes)/ntimes
  matplot(xx_ntimes,t(random_bio),type="l", col="red",xlab='perturbation',ylab="VI",main="VI : Sampled Data")
  matlines(xx_ntimes,t(case_bio),type="l", col="black",xlab='perturbation',ylab="VI",main="VI : Sampled Data")
  plot(ITP.result,main='VI',xrange=c(0,1),xlab='perturbation',ylab="VI")
  lines(xx_ntimes,rep(0.05,21),type="l",col="red")
  
  matrice.quad <- ITP.result$heatmap.matrix[,(p+1):(3*p)]
  ascissa.quad <- seq(abscissa.range[1],abscissa.range[2],length.out=p*2)
  
  image(ascissa.quad,ordinata.grafico,t(matrice.quad[p:1,]),col=colori,ylab='Interval length',main='Adjusted p-value heatmap',xlab='abscissa',zlim=c(0,1),asp=1)
  vertical.image.legend(zlim=c(0,1), col=colori)
  
  dev.off()
  
  
  
  
  ylab='Functional Data'
  lwd=1
  alpha1=0.05
  alpha2=0.01
  xrange=c(0,1)
  main=NULL
  col=c(1,2)
  object=ITP.result
  ylim=range(object$data.eval)
  p <- length(object$pval)
  J <- dim(object$data.eval)[2]
  n <- dim(object$data.eval)[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa.pval = seq(xmin,xmax,len=p)
  Abscissa = seq(xmin,xmax,len=J)
  main.data <- paste(main,': Functional Data')
  main.data <- sub("^ : +", "", main.data)
  colors <- numeric(n)
  colors[which(object$labels==1)] <- "red"
  colors[which(object$labels==2)] <- "blue"
  difference1 <- which(pvalues_ad<alpha1 & pvalues_ad>=alpha2)
  
  #pdf(file_out_plot1)
  setEPS()
  postscript(file_out_plot1)
  
  xx_ntimes=c(0:ntimes)/ntimes
  matplot(xx_ntimes,t(random_bio),type="l", col="white",xlab='perturbation',ylab="VI",main="VI : Sampled Data",cex.lab=2,cex.main=3)
  
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  difference2 <- which(pvalues_ad<alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray80", density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  matlines(xx_ntimes,t(random_bio),type="l", col="blue",xlab='perturbation',ylab="VI",main="VI : Sampled Data",cex.lab=2,cex.main=3)
  matlines(xx_ntimes,t(case_bio),type="l", col="red",xlab='perturbation',ylab="VI",main="VI : Sampled Data")
  
  legend(0.2,0.3,c(as.expression(~VIc[random]),as.expression(~VIc)),lty=c(1,1),lwd=c(2,2),col=c("blue","red"),cex=3)
  dev.off()
  
  ### functional data
  
#  pdf(file_out_plot2)
  setEPS()
  postscript(file_out_plot2)
  
  xx_ntimes=c(0:ntimes)/ntimes
  matplot(xx_ntimes,t(random_bio),type="l", col="white",xlab='perturbation',ylab="VI",main="VI : Data",cex.lab=2,cex.main=3)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  difference2 <- which(pvalues_ad<alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray80", density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  matlines(Abscissa,t(object$data.eval),type='l',main=main.data,ylab=ylab,col=colors,lwd=lwd,ylim=ylim)
  legend(0.45,0.13,c(as.expression(~VIc[random]),as.expression(~VIc)),lty=c(1,1),lwd=c(2,2),col=c("blue","red"),cex=2.3)
  
  dev.off()
  
  
  
  ####copiato dal codice fdatest: plot.ITP2.r
  
  #pdf(file_out_plot3)
  setEPS()
  postscript(file_out_plot3)
  
  main.p <- paste(main,': Adjusted p-values')
  main.p <- sub("^ : +", "", main.p)
  plot(abscissa.pval,pvalues_ad,pch=16,ylim=c(0,1),main=main.p,ylab='p-value',xlab='perturbation',cex.lab=2,cex.main=3)
  
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min.rect <- abscissa.pval[difference1[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray90", density = -2, border = NA)
    }
    #gray90
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  difference2 <- which(pvalues_ad<alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min.rect <- abscissa.pval[difference2[j]] - (abscissa.pval[2] - abscissa.pval[1])/2
      max.rect <- min.rect + (abscissa.pval[2] - abscissa.pval[1])
      rect(min.rect, par("usr")[3], max.rect, par("usr")[4], col = "gray80", density = -2, border = NA)
    }
    #gray80
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL, border = "black")
  }
  for(j in 0:10){
    abline(h=j/10,col='lightgray',lty="dotted")
  }
  points(abscissa.pval,pvalues_ad,pch=16,cex=2)
  lines(xx_ntimes,rep(0.05,21),type="l",col="red")
  
  dev.off()
  
  #######################################
  
}
