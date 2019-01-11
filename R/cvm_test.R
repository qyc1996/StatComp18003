#' @title Cram′er-von Mises test
#' @description  generating Cram′er-von Mises test p-value,which the null hypothesis that the two samples come from the same distribution.
#' @param N the number of permutation
#' @param x the sample to be tested
#' @param y the sample to be tested
#' @return  test p-value
#' @examples
#' \dontrun{
#' set.seed(2)
#' attach(chickwts)
#' x1 <- sort(as.vector(weight[feed == "soybean"]));
#' y1 <- sort(as.vector(weight[feed == "linseed"]));
#' cvm1<-cvm(x1,y1,999);
#' p1 <- mean(cvm1 >= cvm1[1]);
#' p1
#' hist(cvm1, main = "Cram′er-von Mises test", freq = FALSE, xlab = "w (p = 0.385)",breaks = "scott");
#' points(cvm1[1], 0, cex = 1, pch = 16);
#' }
#' @export
cvm<-function(x,y,N){
  reps<-numeric(N);
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
