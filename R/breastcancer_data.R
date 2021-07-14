#' Screened breast cancer data
#'
#' Details of where the data came from
#'
#' @docType data
#'
#' @usage data(breastcancer_screened)
#'
#' @keywords datasets
#'
#' @references Who et. al. (2010) Etc.
#' (\href{https://www.google.com}{Link name here})
#'
#' @source \href{https://www.google.com}{Source}
#'   \code{# importing the original dataset
#'   dat = read_csv("breastcancer.csv")
#'   computing robust correlations and then screening predictors
#'   abscors = lapply(dat[,-1], function(x){abs(paircorxyf((dat$surtime),x))})
#'   x = dat[,names(sort(unlist(abscors),decreasing = T))[1:200]]
#'   y = dat$surtime
#'   write.csv(data.frame(y,x), file = "breastcancer_screened.csv")}
#'
#' @examples
#' data(breastcancer_screened)
#' x = breastcancer_screened[,-1]
#' y = breastcancer_screened$y
#' \donttest{
#' #showing distribution of correlations
#' corrs = paircorf(as.matrix(x))
#' corrplot(corrs, method = "shade",shade.col = NA, tl.col ="black", tl.srt = 45, order = "FPC",tl.cex = 0.01)
#'
#' #showing outlier detection results via DDC
#' fit1 <- DDC(x)
#' cellMap(fit1$remX, fit1$stdResid, columnlabels = c(1,rep(" ",198),200), rowlabels = c(1,rep(" ",76),78), columnangle = 0,
#'         rowtitle = "Observations", columntitle = "Genes", sizetitles = 2,adjustrowlabels = 0.5, adjustcolumnlabels = 0.5)
#'
#'
#' #variable selection results via pariwise estimator and Gaussrank estimator
#' pairfit = covlasso(x,y,cor.method = "pair", scale.method = "qn", lmin = 0.1,adaptive = T)$betahat_opt
#' gaussrankfit = covlasso(x,y,cor.method = "gaussrank", scale.method = "qn", lmin = 0.1, adaptive = T)$betahat_opt
#' }
"breastcancer_screened"
