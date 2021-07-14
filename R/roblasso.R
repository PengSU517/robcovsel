#' Robust adaptive lasso
#'
#' @param x input design matrix
#' @param y input response vector
#' @param cor.method "norm" or "qn"
#' @param scale.method could be "pearson", "pair" and "gaussrank"
#' @param pda.method "nearpd" or FALSE
#' @param lmin min threshold of eigenvalues
#' @param std If TRUE the robust correlation matrix is used,
#'   if FALSE the robust covariance matrix is used.
#' @param adaptive adaptive regularization penalties
#'
#' @return
#' @export
#'
#' @examples
#' #Setting parameters
#' {
#'   n = 100 #sample size
#'   p = 10 #dimension size
#'   e = 0.05 #contamination rate
#'   r = 0.5 #correlation among predictors
#'   gamma = 4 #magnitude of outliers
#'
#'   mu = rep(10,p) # mean vector of predictors
#'
#'   # covariance matrix of the predictor vector
#'   sigma = diag(rep(5^2,p))
#'   for (i in 1:p) {for (j in 1:p) {
#'     if (i !=j)sigma[i,j] = sqrt(sigma[i,i]*sigma[j,j])*r^abs(i-j)}}
#'
#'
#'   beta = c(1,2,3,4,5,rep(0,p-5))# regression parameter vector
#' }
#'
#'
#' #Generating datasets
#' {
#'   xr = MASS::mvrnorm(n,mu,sigma)
#'   error = rnorm(n,0,1)
#'   y = 20+xr%*%beta + error
#'   bi = apply(matrix(0, nrow = n, ncol = p), 2,
#'              function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
#'
#'   outlier = matrix(0,nrow = n,ncol=p)
#'   for (i in 1:n) {
#'     k = sum(bi[i,])
#'     if(k!=0){
#'       label = as.logical(bi[i,])
#'       muc = mu[label]+sigma[label,!label]%*%solve(sigma[!label,!label])%*%(xr[i,!label]-mu[!label])
#'       sigmac = sigma[label,label]-sigma[label,!label]%*%solve(sigma[!label,!label])%*%sigma[!label,label]
#'       outlier[i,label] = muc-gamma*sign(xr[i,label])*diag(sigmac)
#'     }
#'   }
#'   x = xr*(1-bi)+outlier
#'   pairs(data.frame(y = y,x[,1:3]),pch = 19)
#'
#' }
#'
#' paircorf(x) #GK pairwise correlation
#' gaussrankf(x) #gaussrank correlation
#'
#' fit = covlasso(x,y,cor.method = "gaussrank",scale.method = "norm", pda = "nearpd", lmin = 0.1, ada = T)
#' fit
#'
covlasso = function(x, y,
                    cor.method = c("pearson","pair","gaussrank"),
                    scale.method = c("norm","qn"),
                    pda.method = "nearpd",
                    lmin = 0.1,
                    std = TRUE, adaptive = FALSE){

  x = as.matrix(x)
  n = dim(x)[1]
  p = dim(x)[2]
  data = cbind(y,x)

  if(cor.method=="pearson"){cormatrix = cor(data)}
  if(cor.method=="pair"){cormatrix = paircorf(data)}
  if(cor.method=="gaussrank"){cormatrix = gaussrankf(data)}

  if(scale.method=="norm"){
    scale =  apply(data, 2, sd)
  }

  if(scale.method=="qn"){
    scale = apply(data, 2, function(xvec){robustbase::Qn(xvec)})
  }

  ###NearPD calibration
  nearpdf = function(m, method, lmin){

    if(method == "nearpd"){
      rst = eigen(m)
      neweigen = pmax(rst$values, lmin)
      mnew = (rst$vectors)%*%diag(neweigen)%*%t(rst$vectors)
    }

    if(method == F){
      mnew = m
    }
    return(mnew)
  }

  cormatrix = nearpdf(cormatrix, method = pda.method, lmin = lmin)
  covmatrix = (diag(scale))%*%cormatrix%*%(diag(scale))

  if(std){
    Sigma = cormatrix
  } else {
    Sigma = covmatrix
  }

  svdrst = svd(Sigma)
  Sigmasqrt = (svdrst$u)%*%diag(sqrt(svdrst$d))%*%t(svdrst$v)

  response = as.matrix(Sigmasqrt[,1])
  predictor = Sigmasqrt[,-1]


  if(adaptive){
    betastar = MASS::ginv(t(predictor)%*%predictor) %*%t (predictor) %*% response
    weight = as.numeric(abs(1/betastar))
    Predictor = predictor %*% MASS::ginv(diag(weight))
  }else{
    Predictor = predictor
  }

  fit = lars::lars(Predictor,response, type = "lasso",  intercept = F,normalize = F)
  lambda = fit$lambda
  betatilde = t(as.matrix(fit$beta))

  if(adaptive){
    betahat = betatilde/weight
  }else{betahat = betatilde}

  sigmapf = function(beta){(t(response-predictor%*%beta)%*%(response-predictor%*%beta))*ifelse(std,scale[1]^2,1)}
  sigma2hat = apply(betahat, 2, sigmapf)

  if(std){betahat = (scale[1]*(betahat)/scale[-1])}

  penal = apply(betahat, 2, function(betahat){sum(as.logical(betahat))})
  bic = n*log(sigma2hat) + (penal) * log(n)

  betahat = t(betahat)

  label = which.min(bic)
  lambda_opt = lambda[label]
  betahat_opt = betahat[label,]
  sigma2hat_opt = sigma2hat[label]

  list(lambda = lambda, betahat = betahat, sigma2hat = sigma2hat, bic = bic,
       covmatrix = covmatrix, cormatrix = cormatrix,
       lambda_opt = lambda_opt, betahat_opt = betahat_opt, sigma2hat_opt = sigma2hat_opt)
}
