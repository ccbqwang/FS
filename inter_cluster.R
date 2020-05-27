

require(mvtnorm)
require(foreach)
require(doParallel)

GetInitPij <- function(N,J){
  P <- matrix(1/J,N,J)
  for (i in 1:N){
    P[i,] <- rnorm(J, mean=1/J, sd= 0.5 * 1/J )
    P[i,] <- P[i,]/sum(P[i,])
  }
  
  return (P)
}


SAD1.get_mat <- function (par0, times, traits = 1, options = list()) {
  
  par <- par0
  if (class(par0) == "list") 
    par <- unlist(par0)
  t_len <- length(times)
  SAD.1 <- array(0, dim = c(t_len * traits, t_len * traits))
  for (i0 in 1:traits) for (i1 in 1:traits) {
    if (i0 == i1) 
      for (k0 in 1:t_len) for (k1 in k0:t_len) {
        SAD.1[(i0 - 1) * t_len + k0, (i1 - 1) * t_len + 
                k1] <- par[i0 * 2]^2 * par[i0 * 2 - 1]^(k1 - 
                                                             k0)*((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
        SAD.1[(i0 - 1) * t_len + k1, (i1 - 1) * t_len + 
                k0] <- par[i0 * 2]^2 * par[i0 * 2 - 1]^(k1 - 
                                                             k0)*((1-par[i0 * 2 - 1]^(2*k0))/(1-par[i0 * 2 - 1]^2))
      }
  }
  return(SAD.1)
}




StepE <- function(init_par,x,y,times,omiga,N,P,J){
  
  
  covar <- init_par[1:4]
  curvar_par <- matrix(init_par[-c(1:4)],nrow=2,byrow = T)
  
  for( i in 1:N){
    Fi <- rep(0,J)
    for( j in 1:J){
      mu1 <- EI(curvar_par[1,(2*j-1):(2*j)],times[i,])
      mu2 <- EI(curvar_par[2,(2*j-1):(2*j)],times[i,])
      sad <- SAD1.get_mat(covar,times=times[i,],traits=2)
      Fi[j] <- dmvnorm(c(x[i,],y[i,]),c(mu1,mu2),sad)
    }
    
    OmigaF <- c(omiga %*% Fi)
    P[i,] <- (omiga * Fi) / OmigaF
  }
  return(P)
}

