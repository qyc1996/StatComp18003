#' @title  chi-square test statistic
#' @description give two sample vectors,generate list of chi-square test statistic,degree of freedom and p-value.
#' @param x Sample to be tested
#' @param y Sample to be tested
#' @return chi-square test statistic,p-value and degree of freedom.
#' @examples
#' \dontrun{
#' a=31:35;
#' b=c(31,33,35,37,39);
#' m<-cbind(a,b)
#' chisq.test1(a,b)
#' }
#' @export
chisq.test1<- function(x, y){
  #判断输入值是否符合
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
  #计算
  m <- rbind(x, y)
  margin1 <- rowSums(m)
  margin2 <- colSums(m)
  n <- sum(m)
  me <- tcrossprod(margin1, margin2) / n
  x_stat = sum((m - me)^2 / me)#chisq统计量
  dof <- (length(margin1) - 1) * (length(margin2) - 1)#自由度
  p <- pchisq(x_stat, df = dof, lower.tail = FALSE)#p值
  return(list(x_stat = x_stat, df = dof, `p-value` = p))
}
