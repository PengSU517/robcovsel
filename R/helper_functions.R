#' Gaussrank cor
#'
#' @param x input matrix
#'
#' @return
#' @export
#'
#' @examples
#' x = matrix(rnorm(100), ncol = 10)
#' gaussrankf(x)
gaussrankf = function(x){
  n = dim(x)[1]
  xtilde = apply(x,2,function(xvec){qnorm(rank(xvec)/(n + 1))})
  S = cor(xtilde)
  return(S)
}


#' Pairwise correlation using Qn estimator
#'
#' @param x input matrix
#'
#' @return
#' @export
#'
#' @examples
#' x = matrix(rnorm(100), ncol = 10)
#' paircorf(x)
paircorf = function(x){

  n = dim(x)[1]
  p = dim(x)[2]

  sf = function(x) robustbase::Qn(x)

  scale = apply(x, 2, sf)
  stdx = x%*%diag(1/scale)

  I = diag(1,p)
  R = t(as.matrix(rep(1,p)))

  U = stdx%*%(kronecker(I,R) + kronecker(R,I))
  V = stdx%*%(kronecker(I,R) - kronecker(R,I))

  scaleu = apply(U, 2, sf)
  scalev = apply(V, 2, sf)

  cormatrix = matrix((scaleu^2-scalev^2)/(scaleu^2+scalev^2),nrow = p)
  return(cormatrix)

}


#' Pairwise correlation between x and y
#'
#' @param x input vector
#' @param y input vector
#'
#' @return
#' @export
#'
#' @examples
#' x = rnorm(100)
#' y = rnorm(100)
#' paircorxyf(x,y)
paircorxyf = function(x,y){
  sf = function(x) robustbase::Qn(x)
  scalex = sf(x)
  scaley = sf(y)
  u = x/scalex + y/scaley
  v = x/scalex - y/scaley
  cor = ((sf(u))^2 - (sf(v))^2)/((sf(u))^2 + (sf(v))^2)
  return(cor)
}
