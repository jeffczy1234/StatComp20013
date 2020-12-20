#' @title A two dimensional normal Gibbs sampler 
#' @description two dimensional normal Gibbs sampler
#' @param m the number of repeat
#' @param mu1 the mean of the first norm distribution
#' @param mu2 the mean of the first norm distribution
#' @param sigma1 the standard deviation of the first norm distribution
#' @param sigma2 the standard deviation of the first norm distribution
#' @param rho the correlation coefficient
#' @return a Gibbs sampler of the distribution
#' @examples
#' \dontrun{
#' mcmc.norm(1000,1,2,1,4,0.5)
#' }
#' @import stats
#' @export
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
    y.norm[i]=y.go(x.norm[i-1])
  }
  return(cbind(x.norm,y.norm))
}


#' @title using EM algorithm to estimate GMM 
#' @description using EM algorithm to estimate GMM 
#' @param x1 the sample to be solved
#' @return a list of \code{mu} and \code{sigma}
#' @examples
#' \dontrun{
#' mixnorm.em(as.matrix(faithful))
#' }
#' @import mvtnorm
#' @export
mixnorm.em=function(x1){
  n=nrow(x1)
  p=ncol(x1)
  m=p
  Z <- matrix(c(rep(1, n), rep(0, (m - 1) * n)), n, m)
  mu <- matrix(rep(0, m*p), ncol = p, byrow = T)
  for (i in 1:p){
    mu[,i]=mean(x1[,i])*( 0.5+0:(m-1)/(m-1) )
  }
  sig1 <- matrix(rep(0, m*p),ncol = p, byrow = T)
  for (i in 1:p){
    sig1[i,i]=1
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

#' @title A funciton used for illustration.
#' @name Tn62
#' @description a function used in nn
#' @examples
#' \dontrun{
#'Tn62 = function(z, ix, sizes,k) {
#'  n1 = sizes[1] 
#'  n2 = sizes[2]
#'  n = n1 + n2
#'  if(is.vector(z)) {z = data.frame(z,0)}
#'  z = z[ix, ]
#'  NN = nn2(data=z, k=k+1) 
#'  block1 = NN$nn.idx[1:n1,-1]
#'  block2 = NN$nn.idx[(n1+1):n,-1]
#'  i1 = sum(block1 < n1 + .5)
#'  i2 = sum(block2 > n1+.5)
#'  return((i1 + i2) / (k * n))
#'}
#'}
#' @import Ball
#' @import energy
#' @import boot
#' @import RANN
NULL

#' @title A funciton with nn
#' @name eqdist.nn6
#' @description a function used in nn
#' @examples
#' \dontrun{
#'eqdist.nn6 <- function(z,sizes,k){
#'  boot.obj <- boot(data=z,statistic=Tn62,R=999,
#'                   sim = "permutation", sizes = sizes,k=k)
#'  ts = c(boot.obj$t0,boot.obj$t)
#'  p.value =mean(ts>=ts[1])
#'  list(statistic=ts[1],p.value=p.value)
#'}
#'}
#' @import Ball
#' @import energy
#' @import boot
#' @import RANN
NULL

