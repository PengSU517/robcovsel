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



covlasso = function(x, y, cor.method = "gaussrank", scale.method = "qn", center.method = "median", adaptive = TRUE){

  x = as.matrix(x)
  data = cbind(y,x)

  n = dim(x)[1]
  p = dim(x)[2]

  {
    rst = covf(data, cor.method, scale.method, center.method)
    cormatrix = rst$cormatrix
    scale = rst$scale
    center = rst$center
    covmatrix = rst$covmatrix

  }

  Sigma = cormatrix
  eigenrst = eigen(Sigma)
  sqrteigen = pmax(eigenrst$values,0)
  Sigmasqrt = (eigenrst$vectors)%*%diag(sqrteigen)%*%t(eigenrst$vectors)

  response = as.matrix(Sigmasqrt[,1])
  predictor = Sigmasqrt[,-1]

  if(adaptive){
    if(n>2*p){
      betastar = solve(Sigma[-1, -1])%*%Sigma[-1,1]
    }else{
      betastar = solve(Sigma[-1, -1] + diag(rep(0.2, p)))%*%Sigma[-1,1]
    }
    weight = as.numeric(abs(1/betastar))
    Predictor = predictor%*%MASS::ginv(diag(weight))
  }else{
    Predictor = predictor
  }

  fit = lars::lars(Predictor,response, type = "lasso",  intercept = F,normalize = F)
  betatilde = t(as.matrix(fit$beta))
  lambda = fit$lambda

  if(adaptive){
    betahat = betatilde/weight
  }else{betahat = betatilde}

  #sigmapf = function(beta){(t(response-predictor%*%beta)%*%(response-predictor%*%beta))*(scale[1]^2)}
  betahat = (scale[1]*(betahat)/scale[-1])

  sigmapf = function(beta){covmatrix[1,1] - sum(covmatrix[1,-1]*beta)}
  sigma2hat = apply(betahat, 2, sigmapf)

  # betaridge = scale[1]*(solve(Sigma[-1, -1] + diag(rep(0.2, p)))%*%Sigma[-1,1])/scale[-1]
  sigma2_0 = covmatrix[1,1]

  penal = apply(betahat, 2, function(betahat){sum(as.logical(betahat))})
  bic = (sigma2hat/sigma2_0) + (penal) * log(n)/n

  beta0f = function(betahatvec){((center[1] - sum(center[-1]*betahatvec)))}
  beta0hat = apply(betahat, 2, beta0f)

  label = which.min(bic)
  bic_opt = bic[label]
  lambda_opt = lambda[label]####why need to -1
  beta0hat_opt = beta0hat[label]
  betahat_opt = c(const = beta0hat_opt, beta = as.numeric(betahat[,label]))
  sigma2hat_opt = sigma2hat[label]



  list(lambda = lambda, betahat = betahat, beta0hat = beta0hat, sigma2hat = sigma2hat, penal = penal, bic = bic,
       covmatrix = covmatrix, cormatrix = cormatrix, scale = scale, center = center,
       cor.method = cor.method, scale.method = scale.method, adaptive = adaptive,
       lambda_opt = lambda_opt, beta0hat_opt = beta0hat_opt, betahat_opt = betahat_opt, sigma2hat_opt = sigma2hat_opt, bic_opt = bic_opt)
}


#' Title
#'
#' @param x
#' @param y
#' @param cor.method
#' @param scale.method
#' @param center.method
#' @param adaptive
#'
#' @return
#' @export
#'
#' @examples
copulalasso = function(x, y, cor.method = "gaussrank", scale.method = "qn", center.method = "median", adaptive = TRUE){

  x = as.matrix(x)
  data = cbind(y,x)

  n = dim(x)[1]
  p = dim(x)[2]

  {
    rst = covf(data, cor.method, scale.method, center.method)
    datatilde = rst$xtilde
    cormatrix = rst$cormatrix
    scale = rst$scale
    center = rst$center
    covmatrix = rst$covmatrix

  }

  response = datatilde[,1]
  predictor = datatilde[,-1]

  if(adaptive){
    fit0 = glmnet::cv.glmnet(predictor, response, standardize = FALSE, intercept = FALSE, nlambda = 50, alpha = 0)
    betastar = as.numeric(coef(fit0, s = "lambda.1se"))[-1]
    # betastar = solve(cormatrix[-1, -1] + diag(rep(0.2, p)))%*%cormatrix[-1,1]
    weight = as.numeric(abs(1/betastar))
    Predictor = predictor%*%MASS::ginv(diag(weight))
  }else{
    Predictor = predictor
  }

  fit = glmnet::cv.glmnet(Predictor, response, standardize = FALSE, intercept = FALSE, nlambda = 50)
  betatilde = as.numeric(coef(fit, s = "lambda.1se"))[-1]

  if(adaptive){
    betahat = betatilde/weight
  }else{betahat = betatilde}

  betahat = (scale[1]*(betahat)/scale[-1])
  beta0hat = center[1] - sum(center[-1]*betahat)
  betahat = c(beta0hat, betahat)

  list(betahat = betahat, covmatrix = covmatrix, cormatrix = cormatrix, scale = scale, center = center)
}

