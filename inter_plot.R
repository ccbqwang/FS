
Figure1 <- function(ret,ii=c(1,2,3,4),stage,oii){
  
  
  x <- ret$x-0.01
  y <- ret$y-0.01
  nt <- ret$nt-0.1
  sx <- ret$smp_xm-0.01
  sy <- ret$smp_ym-0.01
  snt <- ret$smp_nt-0.1
  pch1 <- c(1,1,1)
  pch2 <- pch1[stage]
  label <- c("A","B","C","D")
  pdf("Figure_1.pdf",width=8,height=7,fonts=c("serif","Palatino"))
  par(fig=c(0,1,0,1),mfrow=c(2,2))
  par(mar=c(2.5,4,0,0),oma=c(2,1.5,2,1))
  k <- 1
  for(i in ii){
    if(max(x[i,])>max(y[i,])){
      plot(NA,NA,xlab=" ",ylab=" ",xlim=c(nt[i,1]-0.1,nt[i,length(nt[i,])]+0.1),ylim=c(min(x[i,1])-0.1,max(x[i,])+0.1),
           mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
    }else{
      plot(NA,NA,xlab=" ",ylab=" ",xlim=c(nt[i,1]-0.1,nt[i,length(nt[i,])]+0.1),ylim=c(min(y[i,1])-0.1,max(y[i,])+0.1),
           mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
    }
    points(nt[i,],x[i,],col="#3A5FCD",pch=pch2[oii[i,]],cex=2)
    lines(snt[i,],sx[i,],lwd=3,col="#3A5FCD")
    points(nt[i,],y[i,],col="#FF4040",pch=pch2[oii[i,]],cex=2)
    lines(snt[i,],sy[i,],lwd=3,col="#FF4040")
    axis(1,nt[i,],rep("",length(nt[i,])),las=1,cex.axis=1.5,tck=-0.03,lwd=0,lwd.ticks=1,mgp=c(2.5,1.5,0),col="grey")
    mtext(label[k],3,cex=1.7,line=-1.9,adj=0.03)
   
    if(k==3){
      mtext("Niche Index",1,cex=1.7,line=3,adj=1.6)
      mtext("Expression on Cell Types",2,cex=1.7,line=3.4,adj=-4.8)
    }
    k <- k + 1
  }
  dev.off()
  
}


Figure11 <- function(ret,ii=c(1,2,3,4),stage,oii){
  
  
  x <- ret$x-0.01
  y <- ret$y-0.01
  nt <- ret$nt-0.1
  sx <- ret$smp_xm-0.01
  sy <- ret$smp_ym-0.01
  snt <- ret$smp_nt-0.1
  pch1 <- c(1,1,1)
  pch2 <- pch1[stage]
  label <- c("A","B","C","D")
  pdf("Figure_1.pdf",width=8,height=7,fonts=c("serif","Palatino"))
  par(fig=c(0,1,0,1),mfrow=c(2,2))
  par(mar=c(2.5,4,0,0),oma=c(2,1.5,2,1))
  k <- 1
  for(i in ii){
    if(max(x[i,])>max(y[i,])){
      plot(NA,NA,xlab=" ",ylab=" ",xlim=c(nt[i,1]-0.1,nt[i,length(nt[i,])]+0.1),ylim=c(min(x[i,1])-0.1,max(x[i,])+0.1),
           mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
    }else{
      plot(NA,NA,xlab=" ",ylab=" ",xlim=c(nt[i,1]-0.1,nt[i,length(nt[i,])]+0.1),ylim=c(min(y[i,1])-0.1,max(y[i,])+0.1),
           mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
    }
    points(nt[i,],x[i,],col="#3A5FCD",pch=pch2[oii[i,]],cex=2)
    lines(snt[i,],sx[i,],lwd=3,col="#3A5FCD")
    points(nt[i,],y[i,],col="#FF4040",pch=pch2[oii[i,]],cex=2)
    lines(snt[i,],sy[i,],lwd=3,col="#FF4040")
    axis(1,nt[i,],rep("",length(nt[i,])),las=1,cex.axis=1.5,tck=-0.03,lwd=0,lwd.ticks=1,mgp=c(2.5,1.5,0),col="grey")
    mtext(label[k],3,cex=1.7,line=-1.9,adj=0.03)
    
    if(k==3){
      mtext("Niche Index",1,cex=1.7,line=3,adj=1.6)
      mtext("Gene Expression in Cell Types",2,cex=1.7,line=3.4,adj=-1.1)
    }
    k <- k + 1
  }
  dev.off()
  
}



Figure2 <- function(ret,ii=c(1,2,3,4),stage,oii){
  
  
  x <- ret$x-0.01
  y <- ret$y-0.01
  nt <- ret$nt-0.1
  
  pch1 <- c(1,1,1)
  pch2 <- pch1[stage]
  label <- c("A","B","C","D")
  pdf("Figure_2.pdf",width=8,height=7,fonts=c("serif","Palatino"))
  par(fig=c(0,1,0,1),mfrow=c(2,2))
  par(mar=c(2.5,4,0,0),oma=c(2,1.5,2,1))
  k <- 1
  for(i in ii){
    
    xf <- EI(para = ret$parax[i,],x=nt[i,])
    yf <- EI(para = ret$paray[i,],x=nt[i,])
    
    xd <- x[i,]-xf
    yd <- y[i,]-yf
    
    plot(NA,NA,xlab=" ",ylab=" ",xlim=c(min(c(xf,yf))-0.1,max(c(xf,yf))+0.1),
         ylim=c(min(c(xd,yd))-0.05,max(c(xd,yd))+0.05),mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
    
    segments(-10,0,10,0,lwd=2,lty=2,col="grey")
    points(xf,xd,col="#3A5FCD",pch=pch2[oii],cex=2)
    points(yf,yd,col="#FF4040",pch=pch2[oii],cex=2)
    
    axis(1,xf,rep("",length(xf)),las=1,cex.axis=1.5,tck=-0.03,lwd=0,lwd.ticks=1,mgp=c(2.5,1.5,0),col="#3A5FCD")
    axis(1,yf,rep("",length(yf)),las=1,cex.axis=1.5,tck=-0.03,lwd=0,lwd.ticks=1,mgp=c(2.5,1.5,0),col="#FF4040")
    mtext(label[k],3,cex=1.7,line=-1.9,adj=0.03)
    
    if(k==3){
      mtext("Power eqation fitting",1,cex=1.7,line=3,adj=3.6)
      mtext("Residuals",2,cex=1.7,line=3.6,adj=1.6)
    }
    k <- k + 1
  }
  dev.off()
  
}





Figure2_c11 <- function(effect){
  
  
  ###1
  
  smt <- effect$smt-0.1
  smx <- effect$smx-0.01
  smy <- effect$smy-0.01
  indf <- effect$indf-0.01
  inds <- effect$inds-0.01
  interf <- effect$interf-0.01
  inters <- effect$inters-0.01
  
  index1 <- which(effect$type==1)
  index2 <- which(effect$type==2)
  index3 <- which(effect$type==3)
  index4 <- which(effect$type==4)
  index5 <- which(effect$type==5)
  
  
  
  
  
  
  pdf("Figure_2c11.pdf",width=6.1,height=12,fonts=c("serif","Palatino"))
  par(fig=c(0,1,0,1),mfrow=c(5,2))
  par(mar=c(2.5,0.5,0,0),oma=c(2.5,7,3,1))
  
  
  #A
  ia1 <- 1;ia2 <- 12
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,2.6),ylim=c(-0.52,1.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index1[ia1],],effect$smx[index1[ia1],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index1[ia1],],effect$indf[index1[ia1],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index1[ia1],],effect$interf[index1[ia1],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index1[ia1],],effect$smy[index1[ia1],],lwd=2,col="red")
  lines(effect$smt[index1[ia1],],effect$inds[index1[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index1[ia1],],effect$inters[index1[ia1],],lwd=2,lty=3,col="red")
  
  text(0.2+2.4*0.15,-0.52+2.12*0.92,effect$gname[index1[ia1]],cex=1.5,font=3,col="black")
  
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("Gene 1",3,cex=1.5,line=0.6)
  mtext("A",3,cex=1.3,line=-1.4,adj=-0.24)
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,2.6),ylim=c(-0.52,1.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index1[ia2],],effect$smx[index1[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index1[ia2],],effect$indf[index1[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index1[ia2],],effect$interf[index1[ia2],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index1[ia2],],effect$smy[index1[ia2],],lwd=2,col="red")
  lines(effect$smt[index1[ia2],],effect$inds[index1[ia2],],lwd=2,lty=2,col="red")
  lines(effect$smt[index1[ia2],],effect$inters[index1[ia2],],lwd=2,lty=3,col="red")
  text(0.2+2.4*0.15,-0.52+2.12*0.92,effect$gname[index1[ia2]],cex=1.5,font=3,col="black")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("Gene 2",3,cex=1.5,line=0.6)
  
  
  #B
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,3.3),ylim=c(-0.52,2.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index2[ia1],],effect$smx[index2[ia1],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index2[ia1],],effect$indf[index2[ia1],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index2[ia1],],effect$interf[index2[ia1],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index2[ia1],],effect$smy[index2[ia1],],lwd=2,col="red")
  lines(effect$smt[index2[ia1],],effect$inds[index2[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index2[ia1],],effect$inters[index2[ia1],],lwd=2,lty=3,col="red")
  text(0.2+3.1*0.2,-0.52+3.12*0.92,effect$gname[index2[ia1]],cex=1.5,font=3,col="black")
  
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("B",3,cex=1.3,line=-1.4,adj=-0.24)
  
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,3.3),ylim=c(-0.52,2.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index2[ia1],],effect$smy[index2[ia1],],lwd=2,col="red")
  lines(effect$smt[index2[ia1],],effect$inds[index2[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index2[ia1],],effect$inters[index2[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index2[ia2],],effect$smx[index2[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index2[ia2],],effect$indf[index2[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index2[ia2],],effect$interf[index2[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.2+3.1*0.2,-0.52+3.12*0.92,effect$gname[index2[ia2]],cex=1.5,font=3,col="black")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  
  
  
  #C
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.6,3.2),ylim=c(-1.12,3.1),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index3[ia1],],effect$smx[index3[ia1],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index3[ia1],],effect$indf[index3[ia1],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index3[ia1],],effect$interf[index3[ia1],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index3[ia1],],effect$smy[index3[ia1],],lwd=2,col="red")
  lines(effect$smt[index3[ia1],],effect$inds[index3[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index3[ia1],],effect$inters[index3[ia1],],lwd=2,lty=3,col="red")
  text(0.6+2.6*0.2,-1.12+4.32*0.92,effect$gname[index3[ia1]],cex=1.5,font=3,col="black")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("C",3,cex=1.3,line=-1.4,adj=-0.24)
  mtext("Expression on Cell Types",2,cex=1.5,line=5)
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.6,3.2),ylim=c(-1.12,3.1),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  
  
  lines(effect$smt[index3[ia2],],effect$smx[index3[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index3[ia2],],effect$indf[index3[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index3[ia2],],effect$interf[index3[ia2],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index3[ia2],],effect$smy[index3[ia2],],lwd=2,col="red")
  lines(effect$smt[index3[ia2],],effect$inds[index3[ia2],],lwd=2,lty=2,col="red")
  lines(effect$smt[index3[ia2],],effect$inters[index3[ia2],],lwd=2,lty=3,col="red")
  text(0.6+2.6*0.2,-1.12+4.32*0.92,effect$gname[index3[ia2]],cex=1.5,font=3,col="black")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  
  
  #D
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.3,2.7),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index4[ia1],],effect$smx[index4[ia1],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index4[ia1],],effect$indf[index4[ia1],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index4[ia1],],effect$interf[index4[ia1],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index4[ia1],],effect$smy[index4[ia1],],lwd=2,col="red")
  lines(effect$smt[index4[ia1],],effect$inds[index4[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index4[ia1],],effect$inters[index4[ia1],],lwd=2,lty=3,col="red")
  
  text(0.3+2.4*0.2,-1.12+3.32*0.92,effect$gname[index4[ia1]],cex=1.5,font=3,col="black")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("D",3,cex=1.3,line=-1.4,adj=-0.24)
  
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.3,2.7),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index4[ia2],],effect$smx[index4[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index4[ia2],],effect$indf[index4[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index4[ia2],],effect$interf[index4[ia2],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index4[ia2],],effect$smy[index4[ia2],],lwd=2,col="red")
  lines(effect$smt[index4[ia2],],effect$inds[index4[ia2],],lwd=2,lty=2,col="red")
  lines(effect$smt[index4[ia2],],effect$inters[index4[ia2],],lwd=2,lty=3,col="red")
  text(0.3+2.4*0.2,-1.12+3.32*0.92,effect$gname[index4[ia2]],cex=1.5,font=3,col="black")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  
  #E
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.4,2.9),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index5[ia1],],effect$smx[index5[ia1],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index5[ia1],],effect$indf[index5[ia1],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index5[ia1],],effect$interf[index5[ia1],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index5[ia1],],effect$smy[index5[ia1],],lwd=2,col="red")
  lines(effect$smt[index5[ia1],],effect$inds[index5[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index5[ia1],],effect$inters[index5[ia1],],lwd=2,lty=3,col="red")
  text(0.3+2.4*0.2,-1.12+3.32*0.92,effect$gname[index5[ia1]],cex=1.5,font=3,col="black")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("E",3,cex=1.3,line=-1.4,adj=-0.24)
  mtext("Niche Index",1,cex=1.5,line=3.2,adj=1.6)
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.4,2.9),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  
  lines(effect$smt[index5[ia2],],effect$smx[index5[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index5[ia2],],effect$indf[index5[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index5[ia2],],effect$interf[index5[ia2],],lwd=2,lty=3,col="#3A5FCD")
  
  lines(effect$smt[index5[ia2],],effect$smy[index5[ia2],],lwd=2,col="red")
  lines(effect$smt[index5[ia2],],effect$inds[index5[ia2],],lwd=2,lty=2,col="red")
  lines(effect$smt[index5[ia2],],effect$inters[index5[ia2],],lwd=2,lty=3,col="red")
  text(0.4+2.5*0.2,-1.12+3.32*0.92,effect$gname[index5[ia2]],cex=1.5,font=3,col="black")
  
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  dev.off()
  
}



Figure2_c1 <- function(effect){
  
  ###1
  
  smt <- effect$smt-0.1
  smx <- effect$smx-0.01
  smy <- effect$smy-0.01
  indf <- effect$indf-0.01
  inds <- effect$inds-0.01
  interf <- effect$interf-0.01
  inters <- effect$inters-0.01
  
  index1 <- which(effect$type==1)
  index2 <- which(effect$type==2)
  index3 <- which(effect$type==3)
  index4 <- which(effect$type==4)
  index5 <- which(effect$type==5)
  
  
  
  
  
  
  pdf("Figure_2c1.pdf",width=6.1,height=12,fonts=c("serif","Palatino"))
  par(fig=c(0,1,0,1),mfrow=c(5,2))
  par(mar=c(2.5,0.5,0,0),oma=c(2.5,7,3,1))
  
  
  #A
  ia1 <- 1;ia2 <- 12
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,2.6),ylim=c(-0.52,1.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index1[ia1],],effect$smx[index1[ia1],],lwd=2,col="red")
  lines(effect$smt[index1[ia1],],effect$indf[index1[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index1[ia1],],effect$interf[index1[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index1[ia2],],effect$smx[index1[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index1[ia2],],effect$indf[index1[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index1[ia2],],effect$interf[index1[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.2+2.4*0.15,-0.52+2.12*0.92,effect$gname[index1[ia1]],cex=1.5,font=3,col="red")
  text(0.2+2.4*0.4,-0.52+2.12*0.92,effect$gname[index1[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("Oocyte",3,cex=1.5,line=0.6)
  mtext("A",3,cex=1.3,line=-1.4,adj=-0.24)
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,2.6),ylim=c(-0.52,1.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index1[ia1],],effect$smy[index1[ia1],],lwd=2,col="red")
  lines(effect$smt[index1[ia1],],effect$inds[index1[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index1[ia1],],effect$inters[index1[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index1[ia2],],effect$smy[index1[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index1[ia2],],effect$inds[index1[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index1[ia2],],effect$inters[index1[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.2+2.4*0.15,-0.52+2.12*0.92,effect$gname[index1[ia1]],cex=1.5,font=3,col="red")
  text(0.2+2.4*0.4,-0.52+2.12*0.92,effect$gname[index1[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("Soma",3,cex=1.5,line=0.6)
  
  
  #B
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,3.3),ylim=c(-0.52,2.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index2[ia1],],effect$smx[index2[ia1],],lwd=2,col="red")
  lines(effect$smt[index2[ia1],],effect$indf[index2[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index2[ia1],],effect$interf[index2[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index2[ia2],],effect$smx[index2[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index2[ia2],],effect$indf[index2[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index2[ia2],],effect$interf[index2[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.2+3.1*0.2,-0.52+3.12*0.92,effect$gname[index2[ia1]],cex=1.5,font=3,col="red")
  text(0.2+3.1*0.55,-0.52+3.12*0.92,effect$gname[index2[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("B",3,cex=1.3,line=-1.4,adj=-0.24)
  
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.2,3.3),ylim=c(-0.52,2.6),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index2[ia1],],effect$smy[index2[ia1],],lwd=2,col="red")
  lines(effect$smt[index2[ia1],],effect$inds[index2[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index2[ia1],],effect$inters[index2[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index2[ia2],],effect$smy[index2[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index2[ia2],],effect$inds[index2[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index2[ia2],],effect$inters[index2[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.2+3.1*0.2,-0.52+3.12*0.92,effect$gname[index2[ia1]],cex=1.5,font=3,col="red")
  text(0.2+3.1*0.55,-0.52+3.12*0.92,effect$gname[index2[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  
  
  
  #C
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.6,3.2),ylim=c(-1.12,3.1),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index3[ia1],],effect$smx[index3[ia1],],lwd=2,col="red")
  lines(effect$smt[index3[ia1],],effect$indf[index3[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index3[ia1],],effect$interf[index3[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index3[ia2],],effect$smx[index3[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index3[ia2],],effect$indf[index3[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index3[ia2],],effect$interf[index3[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.6+2.6*0.2,-1.12+4.32*0.92,effect$gname[index3[ia1]],cex=1.5,font=3,col="red")
  text(0.6+2.6*0.48,-1.12+4.32*0.92,effect$gname[index3[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("C",3,cex=1.3,line=-1.4,adj=-0.24)
  mtext("Expression on Cell Types",2,cex=1.5,line=5)
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.6,3.2),ylim=c(-1.12,3.1),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index3[ia1],],effect$smy[index3[ia1],],lwd=2,col="red")
  lines(effect$smt[index3[ia1],],effect$inds[index3[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index3[ia1],],effect$inters[index3[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index3[ia2],],effect$smy[index3[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index3[ia2],],effect$inds[index3[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index3[ia2],],effect$inters[index3[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.6+2.6*0.2,-1.12+4.32*0.92,effect$gname[index3[ia1]],cex=1.5,font=3,col="red")
  text(0.6+2.6*0.48,-1.12+4.32*0.92,effect$gname[index3[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  
  
  #D
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.3,2.7),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index4[ia1],],effect$smx[index4[ia1],],lwd=2,col="red")
  lines(effect$smt[index4[ia1],],effect$indf[index4[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index4[ia1],],effect$interf[index4[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index4[ia2],],effect$smx[index4[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index4[ia2],],effect$indf[index4[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index4[ia2],],effect$interf[index4[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.3+2.4*0.2,-1.12+3.32*0.92,effect$gname[index4[ia1]],cex=1.5,font=3,col="red")
  text(0.3+2.4*0.38,-1.12+3.32*0.92,effect$gname[index4[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("D",3,cex=1.3,line=-1.4,adj=-0.24)
  
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.3,2.7),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index4[ia1],],effect$smy[index4[ia1],],lwd=2,col="red")
  lines(effect$smt[index4[ia1],],effect$inds[index4[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index4[ia1],],effect$inters[index4[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index4[ia2],],effect$smy[index4[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index4[ia2],],effect$inds[index4[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index4[ia2],],effect$inters[index4[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.7+2.4*0.2,-1.12+3.32*0.92,effect$gname[index4[ia1]],cex=1.5,font=3,col="red")
  text(0.7+2.4*0.38,-1.12+3.32*0.92,effect$gname[index4[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  
  #E
  ia1 <- 2;ia2 <- 4
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.4,2.9),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",cex.axis=1.5,las=1)
  
  lines(effect$smt[index5[ia1],],effect$smx[index5[ia1],],lwd=2,col="red")
  lines(effect$smt[index5[ia1],],effect$indf[index5[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index5[ia1],],effect$interf[index5[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index5[ia2],],effect$smx[index5[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index5[ia2],],effect$indf[index5[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index5[ia2],],effect$interf[index5[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.3+2.4*0.2,-1.12+3.32*0.92,effect$gname[index5[ia1]],cex=1.5,font=3,col="red")
  text(0.3+2.4*0.48,-1.12+3.32*0.92,effect$gname[index5[ia2]],cex=1.55,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  mtext("E",3,cex=1.3,line=-1.4,adj=-0.24)
  mtext("Expression Index",1,cex=1.5,line=3.2,adj=2.6)
  plot(NA,NA,xlab=" ",ylab=" ",xlim=c(0.4,2.9),ylim=c(-1.12,2.2),
       mgp = c(2.7, 1, 0),xaxs="i", yaxs="i",yaxt="n",cex.axis=1.5,las=1)
  
  lines(effect$smt[index5[ia1],],effect$smy[index5[ia1],],lwd=2,col="red")
  lines(effect$smt[index5[ia1],],effect$inds[index5[ia1],],lwd=2,lty=2,col="red")
  lines(effect$smt[index5[ia1],],effect$inters[index5[ia1],],lwd=2,lty=3,col="red")
  
  lines(effect$smt[index5[ia2],],effect$smy[index5[ia2],],lwd=2,col="#3A5FCD")
  lines(effect$smt[index5[ia2],],effect$inds[index5[ia2],],lwd=2,lty=2,col="#3A5FCD")
  lines(effect$smt[index5[ia2],],effect$inters[index5[ia2],],lwd=2,lty=3,col="#3A5FCD")
  text(0.4+2.5*0.2,-1.12+3.32*0.92,effect$gname[index5[ia1]],cex=1.5,font=3,col="red")
  text(0.4+2.5*0.48,-1.12+3.32*0.92,effect$gname[index5[ia2]],cex=1.5,font=3,col="#3A5FCD")
  segments(-100,0,100,0,lwd=1,lty=1,col="black")
  dev.off()
  
}



GO_extract <- function(GO1,thre=0.0001,count=15){
  
  
  n <- length(GO1)
  
  allk <-c()
  for(i in 1:n){
    
    ex1 <- GO1[[i]]
    indd <- (ex1[,6]<thre)+(ex1$Count>count)
    ex11 <- ex1[which(indd==2),3]
    allk <- c(allk,ex11)
    
  }
  allk1 <- unique(allk)
  
  alld <- c()
  for(i in 1:n){
    ex1 <- GO1[[i]]
    alld <- cbind(alld,ex1[match(allk1,ex1[,3]),6])
  }
  
  alld[which(is.na(alld))] <- 1
  rownames(alld) <- allk1
  colnames(alld) <- 1:n
  
  return(alld)
  
}


GO1_extract <- function(GO1,GO2){
  
  
  n <- length(GO1)
  allk <-c()
  for(i in 1:n){
    
    ex1 <- GO1[[i]]
    #indd <- (ex1[,6]<thre)+(ex1$Count>count)
    ex11 <- ex1[which(ex1[,12]==1),4]
    #ex11 <- ex1[which(ii[[i]],3]
    allk <- c(allk,as.character(ex11))
  }
  allk1 <- unique(allk)
  
  alld <- c()
  for(i in 1:n){
    ex1 <- GO2[[i]]
    alld <- cbind(alld,ex1[match(allk1,ex1[,3]),6])
  }
  
  alld[which(is.na(alld))] <- 1
  rownames(alld) <- allk1
  colnames(alld) <- 1:n
  
  return(alld)
  
}








Figure2_c2 <- function(GO1,GO2){
  
  
  #go1 <- GO_extract(GO,thre=0.01,count=5)
  go1 <- GO1_extract(GO1,GO2)
  
  require(pheatmap)
  require(RColorBrewer)
  
  pv <- -log10(go1)
  colnames(pv) <- c("A","B","C","D","E")
  pheatmap(pv,treeheight_row=0, treeheight_col=0,cluster_row = TRUE, cluster_col = FALSE,
           display_numbers = FALSE,colorRampPalette(c("#D6D6D6","#FFEC8B", "#FA8072"))(100),
           angle_col=0,filename = "Figure2_c2.pdf",width=7,height=12,fontsize = 13,fontsize_row=9,
           legend=T,legend_breaks = 0:4, legend_labels = c("0","1", "2", "3", "4"),
           main=expression(paste( "-log"[10], "(", italic(P), ")", sep="" )))
  
}

