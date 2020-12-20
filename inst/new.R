## -----------------------------------------------------------------------------
mcmc.norm=function(m,mu1,mu2,sigma1,sigma2,rho){
  x.go=function(y){
    return( rnorm(1,mu1+rho*sigma1*(y-mu2)/sigma2,(1-rho^2)*sigma1^2) ) 
  }
  y.go=function(x){
    return( rnorm(1,mu2+rho*sigma2*(x-mu1)/sigma1,(1-rho^2)*sigma2^2) ) 
  }
  x.norm=numeric(m)
  y.norm=numeric(m)
  x.norm[1]=mu1
  y.norm[1]=mu2
  for (i in 2:m){
    x.norm[i]=x.go(y.norm[i-1])
    y.norm[i]=y.go(x.norm[i])
  }
  return(cbind(x.norm,y.norm))
}
a.mcmc=mcmc.norm(20000,1,2,1,1,0.5)
cor(a.mcmc[,1],a.mcmc[,2])## the rho is 0.5

## -----------------------------------------------------------------------------
library(mvtnorm)
mixnorm.em=function(x1){
  n=nrow(x1)
  p=ncol(x1)
  m=p
  Z <- matrix(c(rep(1, n), rep(0, (m - 1) * n)), n, m)
  mu <- matrix(rep(0, m*p), ncol = p, byrow = T)
  for (i in 1:p){
    mu[,i]=mean(x1[,i])*( 0.5+0:(m-1)/(m-1) ) ##initial value for mu
  }
  sig1 <- matrix(rep(0, m*p),ncol = p, byrow = T)
  for (i in 1:p){
    sig1[i,i]=1 ##initial value for sigma
  }
  Sig <- matrix(rep(sig1, m),ncol = p, byrow = T)
  pi <- rep(1 / m, m)
  ## EM
  for (t in 1:20) {
    for (k in 1:n) {
      for (l in 1:m) {
        Z[k,l]=pi[l] * dmvnorm(x1[k,], mu[l,], Sig[(p * (l - 1) + 1):(p * l),])
      }
      Z[k,]=Z[k,] / sum(Z[k,])
    }
    pi=colMeans(Z)
    for (i in 1:m) {
      mu[i,] <- t(Z[,i]) %*% x1 / sum(Z[,i])
      sumsig <- Z[1,i] * (x1[1,] - mu[i,]) %*% t(x1[1,] - mu[i,])
      for (k in 2:n) {
        sumsig <- sumsig + Z[k,i] * (x1[k,] - mu[i,]) %*% t(x1[k,] - mu[i,])
      }
      Sig[(p * (i - 1) + 1):(p * i),] <- sumsig / sum(Z[,i])
    }
  }
  out=list(mu,Sig)
  names(out) <- c("mu", "sigma")
  return(out)
}
mixnorm.em(as.matrix(faithful))

