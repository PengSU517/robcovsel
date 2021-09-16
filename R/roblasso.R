#' Robust adaptive lasso
#'
#' @param x input design matrix
#' @param y input response vector
#' @param cor.method "sd" or "qn"
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


covlasso = function(x, y, cor.method = "pair", scale.method = "qn", pda.method = "nearpd",
                    lmin = NULL, std = T, adaptive = T, cormatrix = NULL, scale = NULL){

  x = as.matrix(x)
  data = cbind(y,x)

  n = dim(x)[1]
  p = dim(x)[2]

  if(is.null(cormatrix)){
    rst = covf(data, cor.method, scale.method, pda.method, lmin)
    cormatrix = rst$cormatrix
    scale = rst$scale
    covmatrix = rst$covmatrix
    lmin = rst$lmin
  }else{
    covmatrix = (diag(scale))%*%cormatrix%*%(diag(scale))
  }



  if(std){Sigma = cormatrix
  }else{Sigma = covmatrix}

  svdrst = svd(Sigma)
  Sigmasqrt = (svdrst$u)%*%diag(sqrt(svdrst$d))%*%t(svdrst$v)

  response = as.matrix(Sigmasqrt[,1])
  predictor = Sigmasqrt[,-1]


  if(adaptive){
    betastar = MASS::ginv(t(predictor)%*%predictor)%*%t(predictor)%*%response
    weight = as.numeric(abs(1/betastar))
    Predictor = predictor%*%MASS::ginv(diag(weight))
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
  ebic = n*log(sigma2hat) + (penal) * log(log(n))/n

  betahat = t(betahat)

  label = which.min(bic)
  bic_opt = bic[label]
  lambda_opt = lambda[label-1]
  betahat_opt = betahat[label,]
  sigma2hat_opt = sigma2hat[label]

  list(lambda = lambda, betahat = betahat, sigma2hat = sigma2hat, bic = bic, ebic = ebic,
       covmatrix = covmatrix, cormatrix = cormatrix, lmin = lmin,
       lambda_opt = lambda_opt, betahat_opt = betahat_opt, sigma2hat_opt = sigma2hat_opt, bic_opt = bic_opt)
}
