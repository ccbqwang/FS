
require(DOSE)
require(clusterProfiler)


GO_enrich_type <- function(effect){
  
  
  gname <- effect$gname
  nc <- length(table(effect$type))
  CU_ALL2 <- list()
  for(i in 1:nc){
    Allgenes = bitr(gname, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
    CU_genes <- bitr(gname[which(effect$type==i)], fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
    CU_ALL <- enrichGO(gene = CU_genes$ENTREZID,OrgDb = org.Hs.eg.db,keyType = 'ENTREZID',
                       ont = "ALL",universe=Allgenes$ENTREZID,pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,
                       readable = TRUE)
    CU_ALL1 <- setReadable(CU_ALL, OrgDb = org.Hs.eg.db)
    CU_ALL2[[i]] <- summary(CU_ALL1)
  }
  
  return(CU_ALL2)
}




effect_i <- function(dat,ret,ntl=50){
  
  para <- ret$pm
  nn <- dim(para)[1]
  gname <- rownames(dat$x)
  
  indf <- interf <-inds <- inters <- smx <- smy <- smt <- matrix(NA,nrow=nn,ncol=ntl)
  
  for(i in 1:nn){
    if(is.na(para[i,1])){
      next
    }else{
      nt <- dat$nt[i,]
      x0 <- dat$smp_xm[i,1]
      y0 <- dat$smp_ym[i,1]
      nt1 <- seq(nt[1],nt[length(nt)],length=ntl)
      tmp <- ind.get_mu1(para[i,c(1,2,5,6)],times=nt1,x0=x0,y0=y0)
      indf[i,] <- tmp[1:ntl]; inds[i,] <- tmp[(ntl+1):(2*ntl)]
      smx[i,] <- EI(para=dat$parax[i,],nt1)
      smy[i,] <- EI(para=dat$paray[i,],nt1)
      smt[i,] <- nt1
      interf[i,] <- smx[i,]-indf[i,]
      inters[i,] <- smy[i,]-inds[i,]
    }
  }
  
  
  type <- rep(NA,nn)
  
  for(i in 1:nn){
    if(any(is.nan(interf[i,]))||any(is.nan(inters[i,]))||is.na(para[i,1])){
      next
    }
    if(interf[i,ntl]>0.1&&inters[i,ntl]>0.1){
      type[i] <- 1
    }
    if(interf[i,ntl]>0.1&&abs(inters[i,ntl])<0.1){
      type[i] <- 2
    }
    
    if(abs(interf[i,ntl])<0.1&&inters[i,ntl]>0.1){
      type[i] <- 2
    }
    
    if(interf[i,ntl] < -0.1&&inters[i,ntl] < -0.1){
      type[i] <- 3
    }
    if(interf[i,ntl] < -0.1&&abs(inters[i,ntl]) < 0.1){
      type[i] <- 4
    }
    if(abs(interf[i,ntl]) < 0.1&&inters[i,ntl] < -0.1){
      type[i] <- 4
    }
    if(interf[i,ntl] > 0.1&&inters[i,ntl] < -0.1){
      type[i] <- 5
    }
    if(interf[i,ntl] < -0.1&&inters[i,ntl] > 0.1){
      type[i] <- 5
    }
    if(abs(interf[i,ntl]) < 0.1&&abs(inters[i,ntl]) < 0.1){
      type[i] <- 6
    }
  }
  return(list(indf=indf,inds=inds,interf=interf,inters=inters,type=type,gname=gname,smx=smx,smy=smy,smt=smt))
}




zscore <- function(x){
  
  sdx <- sd(x)
  mx <- mean(x)
  res <- (x-mx)/sdx
  res
  
}


cluster_her <- function(ret,h=6.2){
  
  res <- cbind(ret$x,ret$y)
  res1 <- t(apply(res,1,zscore))
  rownames(res1) <- colnames(res1) <- NULL
  require(pheatmap)
  set.seed(100)
  A <- pheatmap(res1,cluster_row = TRUE, cluster_col = F,color= colorRampPalette(c("navy", "white", "firebrick3"))(500),
                filename="Figure3_cluster.pdf",main = " A B FGC Soma",width=7,height=10,fontsize = 10)
  memb <- cutree(A$tree_row, k = NULL,h=h)
  ret$cluster <- as.numeric(memb)
  return(ret)
}


allo_fit <- function(x,y,nt,ntl=50){
  
  
  nn <- dim(x)[1]
  
  smp_xm <- matrix(NA,nrow=nn,ncol=ntl)
  smp_ym <- matrix(NA,nrow=nn,ncol=ntl)
  smp_nt <- matrix(NA,nrow=nn,ncol=ntl)
  fitparax <- matrix(NA,nrow=nn,ncol=2)
  fitparay <- matrix(NA,nrow=nn,ncol=2)
  for(i in 1:nn){
    nnt <- seq(nt[i,1],nt[i,length(nt[i,])],length=ntl)
    smp_nt[i,] <- nnt
    
    tmp_x <- EI_optim(pheno=x[i,],ES=nt[i,])
    smp_x <- EI(tmp_x$fpar,x=nnt)
    smp_xm[i,] <- smp_x
    fitparax[i,] <- tmp_x$fpar
    
    tmp_y <- EI_optim(pheno=y[i,],ES=nt[i,])
    smp_y <- EI(tmp_y$fpar,x=nnt)
    smp_ym[i,] <- smp_y
    fitparay[i,] <- tmp_y$fpar
    cat("i=",i,tmp_x$LL,tmp_y$LL,"\n")
  }
  return(list(x=x,y=y,nt=nt,smp_xm=smp_xm,smp_ym=smp_ym,smp_nt=smp_nt,parax=fitparax,paray=fitparay))
  
}




dat_index <- function(f,s){
  
  nn <- dim(f)[1]
  fsn <- c()
  fn <- c()
  sn <- c()
  oii <- c()
  for(i in 1:nn){
    fs <- rbind(f[i,],s[i,])
    tmp_fs <- colSums(fs)
    oi <- order(tmp_fs)
    oii <- rbind(oii,oi)
    fsn <- rbind(fsn,as.numeric(tmp_fs[oi]))
    fn <- rbind(fn,as.numeric(f[i,oi]))
    sn <- rbind(sn,as.numeric(s[i,oi]))
  }
  rownames(fsn) <- rownames(fn) <- rownames(sn) <- rownames(f)
  colnames(fsn) <- colnames(fn) <- colnames(sn) <- colnames(f)
  stage <- c(rep(1,4),rep(2,2),rep(3,9))
  return(list(x=fn+0.01,y=sn+0.01,times=fsn+0.1,orderi=oii,stage=stage))
}

ode_fit <- function(dat){
  
  
  nn <- dim(dat$x)[1]
  value <- rep(NA,nn)
  pm <- matrix(NA,nrow=nn,ncol=8)
  
  for(i in 1:nn){
    
    tmp <- ode_fit_ind(x=dat$x[i,],y=dat$y[i,],nt=dat$nt[i,])
    if(is.na(tmp)||class(tmp)=="try-error"){
      cat("Gene=",i," par= ",NA," value= ",NA,"\n")
      next
    }
    pm[i,] <- tmp$par
    value[i] <- tmp$value
    cat("Gene=",i," par= ",tmp$par," value= ",tmp$value,"\n")
  }
  res <- list(pm=pm,value=value)
  return(res)
}



ode_fit_ind <- function(x,y,nt){
  
  nnt <- seq(nt[1],nt[length(nt)],length=20)
  
  tmp_x <- EI_optim(pheno=x,ES=nt)
  smp_x <- EI(tmp_x$fpar,x=nnt)
  
  tmp_y <- EI_optim(pheno=y,ES=nt)
  smp_y <- EI(tmp_y$fpar,x=nnt)
  
  x0 <- smp_x[1]
  y0 <- smp_y[1]
  
  allpara <- c()
  for(i in 1:5){
    para1 <- c(tmp_x$fpar,0,0,tmp_y$fpar,0,0)*runif(8,0.8,1.2)
    o1 <- try(optim(para1,s.mle,s.y=c(smp_x,smp_y),s.t=nnt,x0=x0,y0=y0,method="BFGS",control=list(maxit=2500,trace=F)),TRUE)
    if(class(o1)=="try-error")
      o1 <- try(optim(para1,s.mle,s.y=c(smp_x,smp_y),s.t=nnt,x0=x0,y0=y0,method="Nelder-Mead",control=list(maxit=2500,trace=F)),TRUE)
    if(class(o1)=="try-error")
      next
    allpara <- rbind(allpara,o1$par)
  }
  if(is.null(allpara)){
    return(NA)
  }else{
    o2 <- try(optim(colMeans(allpara),s.mle,s.y=c(smp_x,smp_y),s.t=nnt,x0=x0,y0=y0,method="BFGS",control=list(maxit=2500,trace=F)),TRUE)
    if(class(o2)=="try-error")
      o2 <- try(optim(colMeans(allpara),s.mle,s.y=c(smp_x,smp_y),s.t=nnt,x0=x0,y0=y0,method="Nelder-Mead",control=list(maxit=2500,trace=F)),TRUE)
    return(o2)
  }
}


ind.mle <- function(s.par,s.y,s.t,x0,y0){
  A <- sum((s.y - ind.get_mu1(s.par,s.t,x0,y0))^2 )
  return(A) 
}


### optim function
EI.sum <- function(para,x,y){
  B <- sum((y-EI(para,x))^2)
  return(B)
}


### allometry scale
EI <- function(para,x){
  
  A <- para[1]*x^para[2]
  return(A)
}

#### parameters estimate by BFGS method#######
EI_optim <- function(pheno,ES){
  
  
  X <- ES
  Y <- pheno
  
  plist <-list();LL <- c()
  for(i in 1:10){
    para <- c(runif(2))
    CC <- try(optim(para,EI.sum,x=X,y=Y,method="BFGS",control=list(maxit=1000,trace=F)),TRUE)
    if(class(CC)=="try-error")
      CC <- try(optim(para,EI.sum,x=X,y=Y,method="Nelder-Mead",control=list(maxit=1000,trace=F)),TRUE)
    LL <- c(LL,CC$value);plist[[i]] <- CC$par
  }
  fpar <- plist[[which(LL==min(LL))[1]]]
  
  return(list(X=X,Y=Y,fpar=fpar,LL=min(LL)))
  
}



s.mle <- function(s.par,s.y,s.t,x0,y0){
  A <- sum((s.y - com.get_mu(s.par,s.t,x0,y0))^2 )
  A
}


com.get_mu <- function(par, times, x0,y0)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      b1 = par[2],
      a12 = par[3],
      b12 = par[4],
      a2 = par[5],
      b2 = par[6],
      a21 = par[7],
      b21 = par[8]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f( par0, state0, times );
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- a1*X^b1 + a12*Y^b12
            dY <- a2*Y^b2 + a21*X^b21
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


ind.get_mu1 <- function(par, times, x0,y0)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      a1 = par[1],
      b1 = par[2],
      a2 = par[3],
      b2 = par[4]);
  }
  
  state0 <- c(X=x0, Y=y0);
  y <- COMP.f1( par0, state0, times );
  
  #if (class(y)=="try-error" )
  #  return(rep(NaN, length(times)*2));
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- a1*X^b1
            dY <- a2*Y^b2
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}
