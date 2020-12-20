#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
using namespace Rcpp;

//' @title Generating normal distribution using Rcpp
//' @description Generating normal distribution using Rcpp
//' @param n the number of sample (numeric)
//' @param sigma the standard deviation of the  norm distribution (numeric)
//' @return a normal distribution sample of size \code{n}
//' @examples
//' \dontrun{
//' random(1,100)
//' }
//' @export
// [[Rcpp::export]]
NumericVector random(double sigma, int n) {
  NumericVector runif(n);
  for(int i = 0; i < n; i++) {
    double U, V,Z;
    U = rand()/ (RAND_MAX + 1.0);
    V = rand()/ (RAND_MAX + 1.0);
    Z=sqrt(-2.0 * log(U))* sin(2.0 * 3.141592654 * V);
    runif[i] = sigma*Z;
    
  }
  return runif;
}

#include <Rcpp.h>
#include <stdlib.h>
#include <stdio.h>
using namespace Rcpp;

//' @title Use random walk to get Markov chain  using Rcpp.
//' @description Use random walk to get Markov chain  using Rcpp.
//' @param sig the standard deviation of the  norm distribution (numeric)
//' @param n1 the number of sample (numeric)
//' @return a Markov chain of size \code{n1} and the accept rate
//' @examples
//' \dontrun{
//' chain(1,1000)
//' }
//' @export
// [[Rcpp::export]]

NumericVector chain(double sig, int n1){
  NumericVector rw(n1+1),x1(n1);
  rw[0]=exp(-1)/2;
  x1=random(sig,n1);
  int k1;
  double y1;
  k1=0;
  for(int i = 1; i < n1; i++) {
    y1=random(sig,1)[0]+rw[i-1];
    if (x1[i-1] <= (exp(-abs(y1))/exp(-abs( rw[i-1] ))) )
    {rw[i]=y1;} 
    else {
      rw[i]=rw[i-1];
      k1=k1 + 1;
    }
  }
  rw[n1]=k1;
  return(rw);
}
