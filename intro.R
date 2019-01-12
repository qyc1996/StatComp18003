## ------------------------------------------------------------------------
x<-seq(0,2*pi,length=100)
y1<-sin(x)
y2<-cos(x)
plot(x,y1,type="l",lwd=2,ylab="y")
lines(x,y2,lwd=2,lty=2)
abline(h=0,lwd=2)
abline(v=0,lwd=2)

## ------------------------------------------------------------------------
set.seed(1958);
y<-rnorm(100)
hist(y,freq=FALSE)
res<-density(y);
lines(res,col="red")

## ------------------------------------------------------------------------
library(lattice)
data(quakes); 
mini <- min(quakes$depth) ;
maxi <- max(quakes$depth) ;
int <- ceiling((maxi - mini)/9);
inf <- seq(mini, maxi, int) ;
quakes$depth.cat <- factor(floor(((quakes$depth - mini) / int)), labels=paste(inf, inf + int, sep="-")) ;
xyplot(lat ~ long | depth.cat, data = quakes)

## ------------------------------------------------------------------------
data("pressure");
names(pressure)
lm.summary<-lm(temperature~pressure,data=pressure)
plot(x=pressure$pressure,y=pressure$temperature,xlab="pressure",ylab="temperature")
plot(lm.summary)

## ------------------------------------------------------------------------
cells<-c(1,16,24,68)
rnames<-c("R1","R2")
cnames<-c("C1","C2")
mymatrix<-matrix(cells,nrow=2,ncol=2,byrow=TRUE,dimnames = list(rnames,cnames))
mymatrix

## ------------------------------------------------------------------------
set.seed(2)
x<-c(0,1,2,3,4);
p<-c(0.1,0.2,0.2,0.2,0.3);
cp<-cumsum(p);
m<-10^3;
r<-numeric(m);
r<-x[findInterval(runif(m),cp,left.open = TRUE)+1]#"left.open=TRUE,because in findInterval function ,the intervals are open at right and closed at left.
r
ct<-as.vector(table(r));#calulate frequency of the table "r"
ct/sum(ct);

## ------------------------------------------------------------------------
rnbeta<-function(a,b,N){
  n<-N;
  j<-k<-0;
  y<-numeric(n);
  while(k<n){
    u<-runif(1);
    j<-j+1;
    x<-runif(1) #random variate from g
    if(x^(a-1)*(1-x)^(b-1)>u){
      #we aceept x
 k<-k+1;
 y[k]<-x;}
    
  }
return(y)};
x<-rnbeta(3,2,1000);
p1 = hist(x,plot = T,freq = FALSE,breaks=seq(0,1,0.01),col="yellow",main = "a-r menthod for generating random number");
d=seq(from=min(x),to=max(x),by=0.01) 
lines(x=d,y=dbeta(d,3,2),lty=2,col="red")#add density curve

## ------------------------------------------------------------------------
n<-1000;
r<-4;
beta<-3;
lambda<-rgamma(n,r,beta);
x<-rexp(n,lambda);
hist(x,plot=T,freq=FALSE,breaks=50,col="skyblue",main="Exponential-Gamma mixture")
lines(density(x,bw=1),col='red',lty=2)

## ------------------------------------------------------------------------
x<-seq(0.1,0.9,length = 9)
m<-10000
set.seed(111)
u<-runif(m)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
g<-30*(x[i])^3 *u^2*(1-x[i]*u)^2
cdf[i]<-mean(g) 
}
Phi <- pbeta(x,3,3)
print(round(rbind(x, cdf, Phi), 3))

## ------------------------------------------------------------------------
MC.Phi<-function(x,sigma,R= 10000, antithetic = TRUE) {
u <- runif(R/2)
if (!antithetic)
v<- runif(R/2)
else 
v<-1-u
w<-c(u, v)
cdf <- numeric(length(x))
for (i in 1:length(x)) {
g <-x[i]*w/(sigma)^2*exp((-(x[i]*w)^2/(2*sigma*sigma)))*x[i]
cdf[i] <- mean(g)
}
return(cdf)
}
x <- seq(.1,2.5,length=5);
pRayleigh<-function(x,sigma){
  s<-sigma;
  p<-numeric(length(x));
  intergrand<-function(x){
    x/(s^2)*exp((-x^2/(2*(s^2))))};
  for(i in 1:length(x)){
  p[i]<-integrate(intergrand,0,x[i])$value;
  }
  return(p)
}
Phi<-pRayleigh(x,sigma=2)
set.seed(123)
MC1<- MC.Phi(x,sigma=2,anti=FALSE) #for (X1+X2)/2 which X1,X2 is independent
set.seed(123)
MC2<- MC.Phi(x,sigma=2,anti=TRUE)  #for antithetic variables (X+X')/2
print(round(rbind(x, MC1, MC2, Phi),5))
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 1.95
for (i in 1:m) {
MC1[i] <- MC.Phi(x,2,R = 1000, anti = FALSE)
MC2[i] <- MC.Phi(x,2,R = 1000,anti=TRUE)
}
print(sd(MC1))
print(sd(MC2))
print((var(MC1) - var(MC2))/var(MC1))

## ------------------------------------------------------------------------
x<-seq(1,10,length=100)
y<-1/sqrt(2*pi)*exp(-x^2/2);
z<-exp(-x+1)
plot(x,y,col="red",lty=1,xlim = c(0,10))
lines(x,z,col="green",lty=1,xlim=c(0,10))
y1<-x^2;
z1<-y1*y/z;
plot(x,y1,col="red",lty=1,xlim = c(0,10))
lines(x,z1,col="green",lty=1,xlim=c(0,10))

## ------------------------------------------------------------------------
m <- 10000
set.seed(111)
theta.hat <- se <- numeric(2)
g <- function(x) {
exp(-x^2/2)*x^2/(sqrt(2*pi))*(x >1)
}
x <- rnorm(m) #using f1
fg <- g(x)/dnorm(x)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
set.seed(222)
u <- runif(m) #f2, inverse transform method
x <- 1- log(1 - u)
fg <- g(x) / (exp(-x+1))
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)
print(rbind(theta.hat, se))

## ------------------------------------------------------------------------
n=50
m=1000
G_sam1<-numeric(m)
G_sam2<-numeric(m)
G_sam3<-numeric(m)
set.seed(110)
for(i in 1:m)
{
  x<-sort(rlnorm(n,0,1))
  y<-sort(runif(n,0,1))
  z<-sort(rbinom(n,1,0.5))
  J=2*seq(1:n)-n-1
  G_sam1[i]=(J%*%x)/(n^2*mean(x))
  G_sam2[i]=(J%*%y)/(n^2*mean(y))
  G_sam3[i]=(J%*%z)/(n^2*mean(z))
}
log_norm<-c(mean(G_sam1),median(G_sam1),quantile(G_sam1,seq(0,1,0.1)))
unif_norm<-c(mean(G_sam2),median(G_sam1),quantile(G_sam2,seq(0,1,0.1)))
Bern_norm<-c(mean(G_sam2),median(G_sam1),quantile(G_sam3,seq(0,1,0.1)))
A=rbind(log_norm,unif_norm,Bern_norm)
colnames(A)=c("mean","median",names(quantile(G_sam1,seq(0,1,0.1))))
print(A)
hist(G_sam1,freq=FALSE,col="skyblue",prob=TRUE,main="G density histograms for lognorm")
hist(G_sam2,freq=FALSE,col="skyblue",prob=TRUE,main="G density histograms for  uniform")
hist(G_sam3,freq=FALSE,col="skyblue",prob=TRUE,main="G density histograms for  Bernoulli")

## ------------------------------------------------------------------------
n=20
m=100
G_sam1<-numeric(m)
alpha=0.025
UCL_low<-numeric(1000)
UCL_upp<-numeric(1000)
Mean<-numeric(1000)
for(k in 1:1000)
{ for(i in 1:m) #generate m replitcates
 {
  x<-sort(rlnorm(n,0,1))
  J=2*seq(1:n)-n-1
  G_sam1[i]=(J%*%x)/(n^2*mean(x))
 }

UCL_low[k]<-  mean(G_sam1)-qnorm(1-alpha)*sd(G_sam1)/sqrt(m)
UCL_upp[k]<-  mean(G_sam1)+qnorm(1-alpha)*sd(G_sam1)/sqrt(m)
Mean[k]<-mean(G_sam1)
}
fin_mean<-mean(Mean)
k=sum(UCL_low<fin_mean &UCL_upp>fin_mean)
k/1000

## ------------------------------------------------------------------------
library(MASS)
mu0<-c(1,1)
n<-20
m <- 1000
ptest<-stest<-ktest<-numeric(m)
pho<-seq(0.1,0.9,length=21)
#alternatives
power1 <- numeric(21)
power2 <- numeric(21)
power3 <- numeric(21)
set.seed(122)
for (i in 1:21) {
sigma_i <-matrix(c(1,0.1+0.04*(i-1),0.1+0.04*(i-1),1),nrow = 2)
for(k in 1: m){
x <- mvrnorm(n,mu0,sigma_i)
ptest[k]<- cor.test(x[,1],x[,2],alternative = "greater",method = "pearson")$p.value
stest[k]<-cor.test(x[,1],x[,2],alternative = "greater",method = "spearman")$p.value
ktest[k]<-cor.test(x[,1],x[,2],alternative = "greater",method = "kendall")$p.value
}
power1[i] <- mean(ptest<=0.05)
power2[i] <- mean(stest<=0.05)
power3[i] <- mean(ktest<=0.05)
} 
plot(pho, power1,type="l",col="blue")
lines(pho,power2,type="l",col="red")
lines(pho,power3,type="l",col="yellow")

## ------------------------------------------------------------------------
library(MASS)
m=1000
n=20
mean<-c(1,1)
sigma=matrix(c(1,0,0,1),2,2)
p_pearson<-numeric(m)
p_Spearman<-numeric(m)
p_Kendall<-numeric(m)
for(i in 1:m)
{
  sample=mvrnorm(n,mean,sigma)
  x=sample[,1]
  y=sample[,2]
  test1<-cor.test(x,y,method="pearson")
  test2<-cor.test(x,y,method="kendall")
  test3<-cor.test(x,y,method="spearman")
  p_pearson[i]<-test1$p.value
  p_Spearman[i]<-test2$p.value
  p_Kendall[i]<-test3$p.value
}

ratio=c(mean(p_pearson<0.05),mean(p_Spearman<0.05),mean(p_Kendall<0.05))
names(ratio)<-c("peason","Speqrman","Kendall")
print(ratio)


## ------------------------------------------------------------------------
library(bootstrap);
attach(law);#attach data
n = length( LSAT )
R.hat = cor( LSAT,GPA )
R.jack = numeric ( n )
for (i in 1:n)
  { R.jack[i] = cor( LSAT[-i],GPA[-i] ) } #jackknife
bias.jack = (n-1)*( mean(R.jack) - R.hat )
R.bar = mean(R.jack)
se.jack = sqrt( (n-1)*mean( (R.jack-R.bar)^2 ) )

## ------------------------------------------------------------------------
require(boot); 
attach(aircondit)
x=hours
gaptime.hat = mean(x) #MLE of 1/lambda
B = 5000  #no. boostrap resamples
set.seed(5)
exc75.boot = boot( data=aircondit,statistic=function(x,i){mean(x[i,])}, R=B )

## ------------------------------------------------------------------------
exc75.boot

## ------------------------------------------------------------------------
einf.jack = empinf(exc75.boot, type='jack')#where the call to empinf() is used to generate an empirical influence object using the usual jackknife. 
boot.ci(exc75.boot, type=c('norm','basic','perc','bca'), L=einf.jack)

## ------------------------------------------------------------------------
hist(exc75.boot$t, main='', xlab=expression(1/lambda), prob=T)
points(exc75.boot$t0, 0, pch = 19,col="red")

## ------------------------------------------------------------------------
require(bootstrap); 
attach(scor);#attach data
theta = function(x){
eigen(cov(x))$values[1]/sum(eigen(cov(x))$values)
} #end function
n=length(scor[,1]);
x=as.matrix(scor);
lambda.hat = eigen(cov(x))$values
theta.hat = lambda.hat[1]/sum(lambda.hat)
theta.jack=numeric(n)
for (i in 1:n) {theta.jack[i]=theta( x[-i,] ) }#jackknife
bias.jack=(n-1)*(mean(theta.jack)-theta.hat)
theta.bar=mean(theta.jack)
se.jack = sqrt((n-1)*mean( (theta.jack-theta.bar)^2 ))
print( list(theta.hat=theta(scor), bias=bias.jack, se=se.jack) )
detach(package:bootstrap)

## ------------------------------------------------------------------------
library(DAAG);
attach(ironslag);
n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- 0;
# for n/2-fold cross validation
# fit models on leave-two-out samples
for (k in 1:(n-1)) {
a<-c(k,(k+1))
y <- magnetic[-a]
x <- chemical[-a]
J1 <- lm(y ~ x)
ee1 = magnetic[a]- predict( J1, newdata=data.frame(x=chemical[a]) )
e1[k]<-sum(ee1^2)
J2 <- lm(y ~ x + I(x^2))
ee2 = magnetic[a]- predict( J2, newdata=data.frame(x=chemical[a]) )
e2[k]<-sum(ee2^2)
J3 <- lm(log(y) ~ x)
ee3 = magnetic[a]- exp(predict( J3, newdata=data.frame(x=chemical[a])))
e3[k]<-sum(ee3^2)
J4 <- lm(log(y) ~ log(x))
ee4 = magnetic[a]- exp(predict( J4, newdata=data.frame(x=chemical[a])))
e4[k]<-sum(ee4^2)
}
print( list(mse1=mean(e1/2),mse2=mean(e2/2),mse3=mean(e3/2),mse4=mean(e4/2)) )

## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
cvm<-function(x,y,N){
  reps<-numeric(N);#permutation numbers
  n<-length(x);
  m<-length(y);
  f<-ecdf(x);
  g<-ecdf(y);
  z<-c(x,y);
  t<-numeric(N);
  t0<-m*n/(m+n)^2*(sum((f(x)-g(x))^2)+sum((f(y)-g(y))^2));
  for(i in 1:N){
    k<-sample(length(z),size=n,replace = FALSE);
    x1<-z[k];
    y1<-z[-k];
    f1<-ecdf(x1);
    g1<-ecdf(y1);
    t[i]<-m*n/(m+n)^2*(sum((f1(x)-g1(x))^2)+sum((f1(y)-g1(y))^2))
  } 
  return(c(t0,t));
}
#example 1
set.seed(2)
attach(chickwts)
x1 <- sort(as.vector(weight[feed == "soybean"]));
y1 <- sort(as.vector(weight[feed == "linseed"]));
cvm1<-cvm(x1,y1,999);
p1 <- mean(cvm1 >= cvm1[1]);
p1
hist(cvm1, main = "Cram??er-von Mises test", freq = FALSE, xlab = "w (p = 0.385)",breaks = "scott");
points(cvm1[1], 0, cex = 1, pch = 16);
detach(chickwts)

## ------------------------------------------------------------------------
library(RANN)
library(energy)
library(Ball)
library(boot)
library(ggplot2)
set.seed(123)
m <- 50; k<-3; p<-2;
n1<- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)

Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,mean = 1,sd = 0.6),ncol=p);
  y <- matrix(rnorm(n2*p,mean = 1,sd = 0.8),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ------------------------------------------------------------------------
set.seed(122)
m <- 50; k<-3; p<-2;
n1<- n2 <- 50; R<-999; n <- n1+n2; N = c(n1,n2)

Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- matrix(rnorm(n1*p,mean = 0.3,sd = 0.6),ncol=p);
  y <- matrix(rnorm(n2*p,mean = 0.4,sd = 0.8),ncol=p);
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ------------------------------------------------------------------------
set.seed(124)
m <- 50; k<-3; p<-2; 
n1 <- n2 <- 20; R<-999; n <- n1+n2; N = c(n1,n2)

Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)
for(i in 1:m){
  # t distribution with 1 df (heavy-tailed distribution)
  x <- matrix(rt(n1*p,df = 1),ncol=p); 
  #bimodel distribution (mixture of two normal distributions)
  y <- cbind(rnorm(n2,mean = 0.4),rnorm(n2,mean = 0.5));
  z <- rbind(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ------------------------------------------------------------------------
set.seed(125)
m <- 50; k<-3; p<-2; 
n1 <- 10;n2 <- 100;R<-999; n <- n1+n2; N = c(n1,n2)

Tn <- function(z, ix, sizes,k) {
n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
if(is.vector(z)) z <- data.frame(z,0);
z <- z[ix, ];
NN <- nn2(data=z, k=k+1) # what's the first column?
block1 <- NN$nn.idx[1:n1,-1]
block2 <- NN$nn.idx[(n1+1):n,-1]
i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
(i1 + i2) / (k * n)
}

eqdist.nn <- function(z,sizes,k){
  boot.obj <- boot(data=z,statistic=Tn,R=R,sim = "permutation", sizes = sizes,k=k)
  ts <- c(boot.obj$t0,boot.obj$t)
  p.value <- mean(ts>=ts[1])
  list(statistic=ts[1],p.value=p.value)
}

p.values <- matrix(NA,m,3)
for(i in 1:m){
  x <- c(rnorm(n1,mean = 1,sd = 1)); # n1 = 10
  y <- c(rnorm(n2,mean = 2,sd = 2)); # n2 = 100
  z <- c(x,y)
  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
  p.values[i,3] <- bd.test(x=x,y=y,R=999,seed = i*12)$p.value
}
alpha <- 0.05;
pow <- colMeans(p.values<alpha)
print(pow)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
ggplot(power,aes(methods,pow))+#plot
  geom_col(fill = 'palegreen3')+
  coord_flip()

## ------------------------------------------------------------------------
f903 <- function(x, eta=0, theta=1) {
stopifnot( theta > 0 )
return( 1 / (pi*theta * (1 + ((x-eta)/theta)^2)) )
} #end function
set.seed(903)
m <- 10000
sigma = 1 #try sigma = 1 for proposal scale
x <- numeric(m)
x[1] <- rnorm(1,0,sigma) #initialize with X0~N(0,sigma^2)
k <- 0
u <- runif(m)
for (i in 2:m) {
xt <- x[i-1]
y <- rnorm(1, mean = xt, sd = sigma)
num <- f903(y) * dnorm(xt, mean = y, sd = sigma)
den <- f903(xt) * dnorm(y, mean = xt, sd = sigma)
if (u[i] <= num/den) x[i] <- y else {
x[i] <- xt
k <- k+1 #y is rejected
} #end if/else
} #end for loop
print(k/m) #rejection rate

## ------------------------------------------------------------------------
b0 = 1000 #burn-in
index = (b0+1):m
y1 = x[index] #chain after burn-in

## ------------------------------------------------------------------------
p10 = seq(.1,.9, .1)
round( rbind( quantile(y1, p10), qcauchy(p10) ), 3 )

## ------------------------------------------------------------------------
plot( y1~index, type="l", main="", ylab="x" )

## ------------------------------------------------------------------------
hist( y1, prob=T, main='', xlab='x', ylim=c(0,.35), breaks=50 )
xarg = seq( min(y1), max(y1), 0.1 )
lines( xarg, f903(xarg), ylim=c(0,.35) )

## ------------------------------------------------------------------------
f906 <- function( th,x ) {
if (th<0 || th>1 ) return (0)
(2+th)^x[1] * (1-th)^(x[2]+x[3]) * th^x[4]
} #end function
xdata = c( 125, 18, 20, 34 ) #observed multinom. data
m = 10000
set.seed(906)
th = numeric(m)
th[1] = runif(1) #initialize: sample from prior on theta
k = 0
u = runif(m)
##employ skewed beta proposal density
for (t in 2:m) {
xt = th[t-1]
alph = xt/(1-xt)
y <- rbeta(1, shape1=alph, shape2=1 )
numer = f906( y,xdata ) * dbeta( xt, y/(1-y), 1)
denom = f906( xt,xdata ) * dbeta( y, alph, 1)
if ( u[t] <= numer/denom )
th[t] = y else {
th[t] = th[t-1]
k = k + 1
} #end if/else
} #end for loop

## ------------------------------------------------------------------------
plot( th, type="l", ylim=range(th),
xlab=bquote(theta), ylab="posterior" )

## ------------------------------------------------------------------------
hist( th[2001:m], prob=T, breaks="Scott",
xlab=bquote(theta), ylab='posterior', main="" )

## ------------------------------------------------------------------------
theta.hat = mean( th[2001:m] )
theta.hat

## ------------------------------------------------------------------------
print(c( 0.5 + theta.hat/4, (1 - theta.hat)/4,
(1 - theta.hat)/4, theta.hat/4) ) #p.hat vector


## ------------------------------------------------------------------------
k = c( 4:25, 100, 500, 1000 )
object = function( a, df ){
 a2 = a^2
 arg = sqrt( a2*df/(df + 1 - a2) )
 Sk = pt( q=arg, df=df, lower=F)
 arg = sqrt( a2*(df-1)/(df - a2) )
 Skm1 = pt( q=arg, df=df-1, lower=F)
 return( Sk-Skm1 )
 } #end function
for ( i in 1:length(k) ) {
 print( c( k[i], uniroot(object, lower=1, upper=2, df=k[i])$root ) )
 } #end for loop 

## ------------------------------------------------------------------------
attach(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
#1 for loops
lf_3<- vector("list", length(formulas))
for (i in seq_along(formulas)){
  lf_3[[i]] <- lm(formulas[[i]], data = mtcars)
}
#2 lapply
la_3<-lapply(formulas, function(x) lm(formula = x, data = mtcars))

## ------------------------------------------------------------------------
set.seed(123)
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
# for loops
lf_4<- vector("list", length(bootstraps))
for (i in seq_along(bootstraps)){
  lf_4[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}
# lapply without anonymous function
la_4<- lapply(bootstraps, lm, formula = mpg ~ disp)

## ------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
#For the models in exercise 3:
sapply(la_3, rsq)
sapply(lf_3, rsq)
#For the models in exercise 4:
sapply(la_4,rsq)
sapply(lf_4,rsq)

## ------------------------------------------------------------------------
set.seed(3)
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
# anonymous function:
sapply(trials, function(x) x[["p.value"]])
# without anonymous function:
sapply(trials, "[[", "p.value")

## ------------------------------------------------------------------------
options(warn = -1)
chisq.test1<- function(x, y){
  #?§Ø?????<U+05B5>?<U+01F7>?????
  if (!is.numeric(x)) {
    stop("x must be numeric")}
  if (!is.numeric(y)) {
    stop("y must be numeric")}
  if (length(x) != length(y)) {
    stop("x and y must have the same length")}
  if (length(x) <= 1) {
    stop("length of x must be greater one")}
  if (any(c(x, y) < 0)) {
    stop("all entries of x and y must be greater or equal zero")}
  if (sum(complete.cases(x, y)) != length(x)) {
    stop("there must be no missing values in x and y")}
  if (any(is.null(c(x, y)))) {
    stop("entries of x and y must not be NULL")}
  #????
  m <- rbind(x, y)
  margin1 <- rowSums(m)
  margin2 <- colSums(m)
  n <- sum(m)
  me <- tcrossprod(margin1, margin2) / n
  x_stat = sum((m - me)^2 / me)#chisq<U+0373>??<U+FFFD><U+FFFD>
  dof <- (length(margin1) - 1) * (length(margin2) - 1)#???<U+0276>? 
  p <- pchisq(x_stat, df = dof, lower.tail = FALSE)#p<U+05B5>
  return(list(x_stat = x_stat, df = dof, `p-value` = p))
}
#check
a=31:35;
b=c(31,33,35,37,39);
m<-cbind(a,b)
chisq.test1(a,b)
chisq.test(m)
microbenchmark::microbenchmark(
  chisq.test(m),
  chisq.test1(a, b))

## ------------------------------------------------------------------------
table2 <- function(x, y) {
  x_val <- unique(x)
  y_val <- unique(y)
  mat <- matrix(0L, length(x_val), length(y_val))
  for (i in seq_along(x)) {
    mat[which(x_val == x[[i]]), which(y_val == y[[i]])] <-
      mat[which(x_val == x[[i]]),  which(y_val == y[[i]])] + 1L
  }
  dimnames <- list(x_val, y_val)
  names(dimnames) <- as.character(as.list(match.call())[-1])  # R has names for dimnames... :/
  tab <- array(mat, dim = dim(mat), dimnames = dimnames)
  class(tab) <- "table"
  tab
}
a <- c(1, 2, 3, 1, 2, 3)
b <- c(2, 3, 4, 2, 3, 4)
identical(table(a, b), table2(a, b))
microbenchmark::microbenchmark(table(a, b), table2(a, b))
c<-c(1,2,3,4,5)
identical(table(c, c), table2(c, c))
microbenchmark::microbenchmark(table(c, c), table2(c, c))

