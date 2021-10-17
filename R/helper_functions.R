
#' Correlation estimators
#'
#' @param x input matrix
#' @param cor.method "sd"  or "qn"
#' @param scale.method could be "pearson", "pair" and "gaussrank"
#' @param pda.method "nearpd" or FALSE
#' @param lmin min threshold of eigenvalues
#'
#'
#' @return
#' @export
#'
#' @examples
#' x = matrix(rnorm(100), ncol = 10)
#' paircorf(x)
#'
covf = function(x, cor.method, scale.method, pda.method, lmin = NULL){

  n = dim(x)[1]
  p = dim(x)[2]

  if(scale.method == "sd"){sf = function(x) sd(x)}
  if(scale.method == "mad"){sf = function(x) mad(x)}
  if(scale.method == "qn"){sf = function(x) robustbase::Qn(x)}

  scale = apply(x,2,sf)


  if(cor.method=="pearson"){
    cormatrix = cor(x)
  }

  if(cor.method=="gaussrank"){
    xtilde = apply(x,2,function(xvec){qnorm(rank(xvec)/(n + 1))})
    cormatrix = cor(xtilde)
  }

  if(cor.method== "pair"){

    stdx = x%*%diag(1/scale)
    I = diag(1,p)
    R = t(as.matrix(rep(1,p)))

    U = stdx%*%(kronecker(I,R) + kronecker(R,I))
    V = stdx%*%(kronecker(I,R) - kronecker(R,I))

    scaleu = apply(U, 2, sf)
    scalev = apply(V, 2, sf)

    cormatrix = matrix((scaleu^2-scalev^2)/(scaleu^2+scalev^2),nrow = p)
  }

  if(pda.method!=F){
    nearpdrst = nearpdf(cormatrix, pda.method, lmin)
    cormatrix = nearpdrst$cormatrix.pd
    lmin = nearpdrst$lmin}
  covmatrix = (diag(scale))%*%cormatrix%*%(diag(scale))
  return(list(covmatrix = covmatrix, cormatrix = cormatrix, scale = scale,lmin = lmin))
}



#' PD calibration
#'
#' @param cormatrix input cormatrix
#' @param pda.method "nearpd" or FALSE
#' @param lmin min threshold of eigenvalues
#'
#'
#' @return
#' @export
#'
#' @examples
#' x = matrix(rnorm(100), ncol = 10)
#' nearpdf(cor(x), pda.method = "nearpd", lmin = 0.1)

nearpdf = function(cormatrix, pda.method, lmin){

  if(is.null(lmin)){lmin = lminsel(cormatrix)}

  if(pda.method == "shrinkage"){
    rst = eigen(cormatrix, EISPACK = T)
    delta = min(max((lmin - min(rst$values))/(1- min(rst$values)), 0),1)
    neweigen = delta + (1-delta)*rst$values
    cormatrix.pd = (rst$vectors)%*%diag(neweigen)%*%t(rst$vectors)
  }

  if(pda.method == "nearpd"){
    rst = eigen(cormatrix, EISPACK = T)
    neweigen = pmax(rst$values, lmin)
    cormatrix.pd = (rst$vectors)%*%diag(neweigen)%*%t(rst$vectors)
  }
  return(list(cormatrix.pd = cormatrix.pd, lmin = lmin)  )
}



#' Selection of the eigen value threshold
#'
#' @param cormatrix input cormatrix
#' @return
#' @export
#'

lminsel = function(cormatrix){

  rst = eigen(cormatrix, EISPACK = T)
  lossf = function(lmin, multi = 1){
    neweigen = pmax(rst$values, lmin)
    cormatrixnew = (rst$vectors)%*%diag(neweigen)%*%t(rst$vectors)
    return(matrixcalc::frobenius.norm(cormatrix-cormatrixnew) - multi*log10(lmin))
  }
  lmin = optimize(interval = c(0,2), f =  lossf)$minimum
  return(lmin)

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

#' generating predictors and response randomly using default settings
#'
#' @param n sample size
#' @param p dimension
#' @param e contamination rate
#' @param r correlation coefficent
#' @param gamma magnitude of outliers
#' @param beta
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat = genevar()
#' x = dat$x
#' y = dat$y
#'
genevar = function(n = 100, p = 20, e = 0.05, r = 0.5, beta = c(1,2,1,2,1,rep(0,p-5)), gamma = 10){
  {
    mu = rep(0,p)
    sigma = diag(rep(1^2,p))
    for (i in 1:p) {for (j in 1:p) {
      if (i !=j)sigma[i,j] = sqrt(sigma[i,i]*sigma[j,j])*r^abs(i-j)}}
  }

  {
    xr = mvtnorm::rmvnorm(n = n, mean = mu, sigma = sigma)
    error = rnorm(n,0,0.5)
    y = xr%*%beta + error
    bi = apply(matrix(0, nrow = n, ncol = p), 2,
               function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
    outl = rnorm(n = n*p, mean = gamma, sd = 1)
    rsign = sample(c(-1,1), size = n*p, replace = T)
    outlier = matrix(outl*rsign, nrow = n, ncol=p)
    x = xr*(1-bi)+outlier*bi
  }
  return(list(x = x, y = y, beta = beta))


}







