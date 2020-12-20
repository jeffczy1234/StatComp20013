## -----------------------------------------------------------------------------
u11 <- runif(500)
x11=2/sqrt(1-u11)##分布函数的反函数
hist(x11, prob = TRUE, breaks=seq(0,ceiling(max(x11)),1),main = expression(f(x)== 8*x^(-3) ) ) 
y11 <- seq(0,ceiling(max(x11)),1)
lines(y11, 8*y11^(-3),col = "red")

## -----------------------------------------------------------------------------
y12=c()
for(i in 1:5000){
  x12=runif(3,-1,1)
  if( ( abs(x12[3])>=abs(x12[1]) )&&( abs(x12[3])>=abs(x12[2]) ) )   {y12=c(y12,x12[2])}
  else   {y12=c(y12,x12[3])}
}##生成密度函数的样本
hist(y12, prob = TRUE, breaks=seq(-1,1,0.1),main = expression(f(x)== 3*(1-x^(2))/4 ) ) 
z12 <- seq(-1,1,0.1)
lines(z12, 3/4*(1-z12^(2)),col = "red")

## -----------------------------------------------------------------------------
u14 <- runif(1000)
x14=2/(1-u14)^(1/4)-2##分布函数的反函数
hist(x14, prob = TRUE, breaks=seq(0,ceiling(max(x14)),0.1),main = expression(f(x)== 64*(2+x)^(-5) ) ) 
y14 <- seq(0,ceiling(max(x11)),0.1)
lines(y14, 64*(2+y14)^(-5),col = "red")

## -----------------------------------------------------------------------------
a21=runif(10000, min=0, max=1/3*(pi))
b21=mean(sin(a21)) *1/(1/3*(pi))

## -----------------------------------------------------------------------------
g21<-function(x) {sin(x)}
c21=integrate(g21,0,1/3*(pi))$value

## -----------------------------------------------------------------------------
c(abs(b21-c21),abs(b21-c21)/b21)

## -----------------------------------------------------------------------------
#对偶法
u221=runif(5000)
a221=sum(exp(u221)+exp(1-u221) )/10000
##mc法
u222=runif(10000)
a222=mean(exp(u222))
##积分计算
g22<-function(x) {exp(x)}
e22=integrate(g22,0,1)$value
c(a221,a222,e22)

## -----------------------------------------------------------------------------
b221=var( (exp(u221)+exp(1-u221))/2 )
b222=var(exp(u222))
c(b222-b221,b221/b222)

## -----------------------------------------------------------------------------
##定义拒绝域的生成函数b1
f411<-function(m){
mean411=mean(m)
a411=mean( (m-mean411)^3 )
a412=mean( (m-mean411)^2 )
return( a411/a412^(1.5) )
}
##五个不同的取样次数
n411=c(100,200,500,1000,2000)
cv411 <- qnorm(.975, 0, sqrt(6*(n411-2) / ((n411+1)*(n411+3))))
f412<-function(beta){
k411=numeric(length(n411))
for (i in 1:length(n411)){
  k412=numeric(10000)
  for(j in 1:10000){
    b411=rbeta(n411[i],beta,beta)
    k412[j]=as.integer( abs(f411(b411)) >= cv411[i] )
  }
  k411[i]=mean(k412)
}
return(k411)
}
###选择了五个不同参数
cbind(f412(1),f412(5),f412(10),f412(100),f412(1000))


## -----------------------------------------------------------------------------
##讲随机数换成t分布
f413<-function(beta){
  k411=numeric(length(n411))
  for (i in 1:length(n411)){
    k412=numeric(10000)
    for(j in 1:10000){
      b411=rt(n411[i],beta)
      k412[j]=as.integer( abs(f411(b411)) >= cv411[i] )
    }
    k411[i]=mean(k412)
  }
  return(k411)
}
###选择了五个不同参数
cbind(f413(1),f413(5),f413(10),f413(100),f413(1000))

## -----------------------------------------------------------------------------
##定义Count Five函数
f421=function(x421, y421) {
  X421=x421 - mean(x421)
  Y421=y421 - mean(y421)
  a421<- sum(X421 > max(Y421)) + sum(X421 < min(Y421))
  b421<- sum(Y421 > max(X421)) + sum(Y421 < min(X421))
  return(as.integer(max(c(a421, b421)) > 5))
}
f422=function(x422){
  power1<-replicate(1000, expr={
    x4 <- rnorm(x422, 0, 1)
    y4 <- rnorm(x422, 0, 1.5)
    a411=var.test(x4, y4, ratio = 1,
                  alternative = c("two.sided", "less", "greater"),
                  conf.level = 0.95)$p.value
    c(a411,f421(x4, y4))
  })
  c(mean(power1[1,]),mean(power1[2,]))
}
##第一行是F检验的，第二行是five函数
cbind(f422(6),f422(8),f422(10),f422(20),f422(40),f422(60),f422(80),f422(100),f422(120),f422(140))

## -----------------------------------------------------------------------------
##先重复6.8
##定义好d的求值函数
f431<-function(m){
  mean431=mean(m)
  sigma431=sum( (m-mean431)^2 )/length(m)
  a411=numeric(0)
  a431=t(t( (m-mean431)^3 ))%*% t( (m-mean431)^3 )
  a411=sum(a431)/(length(m)^2*sigma431^3)
  return( a411 )
}
##给定取样次数
n431=c(50,100,200,400,700,1000)
b431=(1/n431)*6*qchisq(0.95,1)
f432=function(beta){
  k431=numeric(length(n431))
  for (i in 1:length(n431)){
    k432=numeric(1000)
    for(j in 1:1000){
      c431=rbeta(n431[i],beta,beta)
      k432[j]=as.integer( abs(f431(c431)) >= b431[i] )
    }
    k431[i]=mean(k432)
  }
  return(k431)
}

cbind(f432(1),f432(5),f432(10),f432(100),f432(1000))

## -----------------------------------------------------------------------------
##复现6.10
n433=30
m433=2500
epsilon433=seq(0, 1, 0.05)
k433=numeric(length(epsilon433))
##使用b31为定值参数
for (j in 1:length(epsilon433)) { 
  d433=epsilon433[j]
  a433<-numeric(m433)
  for (i in 1:m433) { 
    sigma4=sample(c(1, 10), replace = TRUE,size = n433, prob = c(1-d433, d433))
    x=rnorm(n433, 0, sigma4)
    a433[i]=as.integer( abs( f431(x) ) >= (6*qchisq(0.95,1)/n433) )
  }
  k433[j]=mean(a433)
}
#plot power vs epsilon
plot(epsilon433, k433, type = "b",
     xlab = bquote(epsilon433), ylim = c(0,1))
abline(h = 0.05, lty = 3)

## -----------------------------------------------------------------------------
##两个不同的样本量
x61=rnorm(20, 1, sd = 2)
y61=rnorm(15, 1, sd = 2)
##使用count5比较
count5test6 <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(as.integer(max(c(outx, outy)) > 5))
}
R61 = 999;
z61 = c(x61, y61)
K61 = 1:35
set.seed(11111)
rep61 = numeric(R61)
t610 = count5test6(x61, y61)
for (i in 1:R61) {
  k61 =sample(K61, size = length(x61), replace = FALSE)
  x61 = z61[k61]
  y61 = z61[-k61]
  rep61[i] = count5test6(x61, y61)
}
p61=mean(c(t610,rep61))
round(c(p61),3)


## -----------------------------------------------------------------------------
library("Ball")
library("energy")
library(boot)
library(RANN)
Tn62 = function(z, ix, sizes,k) {
  n1 = sizes[1] 
  n2 = sizes[2]
  n = n1 + n2
  if(is.vector(z)) {z = data.frame(z,0)}
  z = z[ix, ]
  NN = nn2(data=z, k=k+1) 
  block1 = NN$nn.idx[1:n1,-1]
  block2 = NN$nn.idx[(n1+1):n,-1]
  i1 = sum(block1 < n1 + .5)
  i2 = sum(block2 > n1+.5)
  return((i1 + i2) / (k * n))
}
eqdist.nn6 <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn62,R=999,
                   sim = "permutation", sizes = sizes,k=k)
  ts = c(boot.obj$t0,boot.obj$t)
  p.value =mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}
##下面分别对不同情况计算p值
##情况一
m6 = 300##迭代次数
set.seed(10001)
n61 = n62 = 50##样本数
k6=3; n6 <- n61+n62; N6 = c(n61,n62)
p6=1#矩阵列数
##方差和均值
sd61 =1
sd62=2
mu61=mu62=0
R6=999##随机次数
p.values61 <- matrix(NA,m6,3)
for(i in 1:m6){
  x6 = matrix(rnorm(n61,mu61,sd61),ncol=p6);
  y6 = cbind(rnorm(n62,mu62,sd62));
  z6 = rbind(x6,y6)
  p.values61[i,1] = eqdist.nn6(z6,N6,k6)$p.value
  p.values61[i,2] = eqdist.etest(z6,sizes=N6,R=R6)$p.value
  p.values61[i,3] = bd.test(x=x6,y=y6,num.permutations=999,seed=i*12345)$p.value
}
alpha61 = 0.05;
pow61 = colMeans(p.values61<alpha61)
round(pow61,5)


## -----------------------------------------------------------------------------
##情况二
m6 = 300##迭代次数
R6=999##随机次数
set.seed(10002)
n61 = n62 = 50##样本数
n6 = n61+n62; N6 = c(n61,n62)
##不同的方差和均值
sd61 =1;sd62=1.5;mu61=1;mu62=0
p.values62 <- matrix(NA,m6,3)
for(i in 1:m6){
  x6 = matrix(rnorm(n61,mu61,sd61),ncol=p6)
  y6 = cbind(rnorm(n62,mu62,sd62));
  z6 = rbind(x6,y6)
  p.values62[i,1] = eqdist.nn6(z6,N6,k6)$p.value
  p.values62[i,2] = eqdist.etest(z6,sizes=N6,R=999)$p.value
  p.values62[i,3] = bd.test(x=x6,y=y6,num.permutations=999,seed=i*12345)$p.value
}
alpha61 = 0.05;
pow62 = colMeans(p.values62<alpha61)
print(pow62)

## -----------------------------------------------------------------------------
##情况三
m6 = 300##迭代次数
R6=999##随机次数
set.seed(11113)
n61 = n62 = 50##样本数
n6 = n61+n62; N6 = c(n61,n62)
p.values63 <- matrix(NA,m6,3)
for(i in 1:m6){
  ##先考虑t分布
  x6 = matrix(rt(n61,1),ncol=p6);
  ##将混合正态的均值也设成0
  y6 = cbind(0.3*rnorm(n62,0,1)+0.7*rnorm(n62,0,2))
  z6 = rbind(x6,y6)
  p.values63[i,1] = eqdist.nn6(z6,N6,k6)$p.value
  p.values63[i,2] = eqdist.etest(z6,sizes=N6,R=R6)$p.value
  p.values63[i,3] = bd.test(x=x6,y=y6,num.permutations=999,seed=i*12345)$p.value
}
alpha61 = 0.05;
pow63 = colMeans(p.values63<alpha61)
print(pow63)

## -----------------------------------------------------------------------------
##情况四
m6 = 300##迭代次数
R6=999##随机次数
set.seed(10002)
n61=50;n62=30##不同的样本数
n6 = n61+n62; N6 = c(n61,n62)
##不同的方差和相同的均值
sd61 =1;sd62=1;mu61=1;mu62=0
p.values64<- matrix(NA,m6,3)
for(i in 1:m6){
  x6 = matrix(rnorm(n61,mu61,sd61),ncol=p6)
  y6 = cbind(rnorm(n62,mu62,sd62));
  z6 = rbind(x6,y6)
  p.values64[i,1] = eqdist.nn6(z6,N6,k6)$p.value
  p.values64[i,2] = eqdist.etest(z6,sizes=N6,R=R6)$p.value
  p.values64[i,3] = bd.test(x=x6,y=y6,num.permutations=999,seed=i*12345)$p.value
}
alpha61 = 0.05;
pow64 = colMeans(p.values64<alpha61)
print(pow64)
##不同的方差和相同的均值
sd61 =1;sd62=2;mu61=0;mu62=0
p.values65<- matrix(NA,m6,3)
for(i in 1:m6){
  x6 = matrix(rnorm(n61,mu61,sd61),ncol=p6)
  y6 = cbind(rnorm(n62,mu62,sd62));
  z6 = rbind(x6,y6)
  p.values65[i,1] = eqdist.nn6(z6,N6,k6)$p.value
  p.values65[i,2] = eqdist.etest(z6,sizes=N6,R=R6)$p.value
  p.values65[i,3] = bd.test(x=x6,y=y6,num.permutations=999,seed=i*12345)$p.value
}
alpha61 = 0.05;
pow65 = colMeans(p.values65<alpha61)
print(pow65)

## ----out.width ='\\textwidth'-------------------------------------------------
library(Rcpp)
##如下，random为生成正态随机数的cpp函数，chain为生成随机游走的函数
sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
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
  
//[[Rcpp::export]]  
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
')
##由于c语言画图的可视化比较复杂，采用R语言本身的画图来比较链的区别
N91=1000
sigma91=c(0.1,0.7, 1,4)
rw911=chain(sigma91[1], N91)
rw912=chain(sigma91[2], N91)
rw913=chain(sigma91[3], N91)
rw914=chain(sigma91[4], N91)
##对比拒绝的比例并画图比较数据
c(rw911[(N91+1)], rw912[(N91+1)], rw913[(N91+1)], rw914[(N91+1)])/N91


## -----------------------------------------------------------------------------
library(Rcpp)
qqplot(random(1,100),rnorm(100,0,1))
qqplot(random(1,1000),rnorm(1000,0,1))
qqplot(random(1,10000),rnorm(10000,0,1))


## -----------------------------------------------------------------------------
library(microbenchmark)
x93=1000
a93=microbenchmark(rnorm=runif(x93,0,2),randomcpp=random(2,x93))
print(summary(a93)[,c(1,3,5,6)])

