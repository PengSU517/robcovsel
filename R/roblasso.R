#' Robust adaptive lasso based on robust correlation estimates
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
#' @param cormatrix you could also use a correlation matrix as input
#' @param scale put in scales if you use cormatrix
#'
#' @return betahat_opt, the optimal estimation of beta obtained from this algorithm
#' @return lambda_opt is the optimal tuning parameter and sigma_opt is the optimal estimation of sigma.
#' @return The output also includes the estimated correlation matrix, the estimated covariance matrix and et cetera from the covf function.
#' @export
#' @examples
#' dat = genevar()
#' y = dat$y
#' x = dat$x
#' fit = covlasso(x,y)
#' fit$betahat_opt
#'



covlasso = function(x, y, cor.method = "gaussrank", scale.method = "qn", center.method = "median", pda.method = NULL,
                    lmin = NULL, std = TRUE, adaptive = TRUE){

  x = as.matrix(x)
  data = cbind(y,x)

  n = dim(x)[1]
  p = dim(x)[2]

  if(is.null(pda.method)){
    if((cor.method == "pair")|(p>(0.5*n))){pda.method = "nearpd"}else{
      pda.method = F
    }
  }

  {
    rst = covf(data, cor.method, scale.method, center.method, pda.method, lmin)
    cormatrix = rst$cormatrix
    scale = rst$scale
    center = rst$center
    covmatrix = rst$covmatrix
    lmin = rst$lmin
  }



  if(std){Sigma = cormatrix
  }else{Sigma = covmatrix}

  eigenrst = eigen(Sigma)
  Sigmasqrt = (eigenrst$vectors)%*%diag(sqrt(eigenrst$values))%*%t(eigenrst$vectors)

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

  #betahat = t(betahat)
  beta0f = function(betahatvec){((center[1] - sum(center[-1]*betahatvec)))}
  #*ifelse(std,scale[1]^2,1)
  beta0hat = apply(betahat, 2, beta0f)

  label = which.min(bic)
  bic_opt = bic[label]
  lambda_opt = lambda[label]####why need to -1
  beta0hat_opt = beta0hat[label]
  betahat_opt = c(const = beta0hat_opt, beta = as.numeric(betahat[,label]))
  sigma2hat_opt = sigma2hat[label]



  list(lambda = lambda, betahat = betahat, beta0hat = beta0hat, sigma2hat = sigma2hat, penal = penal, bic = bic,
       covmatrix = covmatrix, cormatrix = cormatrix, scale = scale, center = center,
       cor.method = cor.method, scale.method = scale.method, pda.method = pda.method,
       lmin = lmin, std = std, adaptive = adaptive,
       lambda_opt = lambda_opt,
       betahat_opt = betahat_opt, sigma2hat_opt = sigma2hat_opt, bic_opt = bic_opt)
}
