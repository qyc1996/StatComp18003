#' @title random numbers which are Rayleighed distribution
#' @description the function for generating distribution function random numbers
#' @param R the number of replication samples
#' @param x the initiation sequence
#' @param sigma the parameter of Rayleigh density
#' @param antithetic the method which generates the cdf.default is TRUE
#' @return a random vector that the number is length(x)
#' @examples
#' \dontrun{
#' x <- seq(.1,2.5,length=5);
#' set.seed(123)
#' MC1<- MC.Phi(x,sigma=2,anti=FALSE) #for (X1+X2)/2 which X1,X2 is independent
#' MC2<- MC.Phi(x,sigma=2,anti=TRUE)  #for antithetic variables (X+X')/2
#' print(round(rbind(x, MC1, MC2, Phi),5))
#' }
#' @export
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
