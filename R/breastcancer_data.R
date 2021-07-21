#' Screened breast cancer data
#'
#' Details of where the data came from
#'
#' @docType data
#'
#' @usage data(breastcancer_screened)
#'
#' @keywords datasets; Breast Cancer
#'
#'
#' @source van 't Veer, L. J., Dai, H., van de Vijver, M. J., He, Y. D., Hart, A. A. M., Mao, M., Peterse,347H. L., van der Kooy, K., Marton, M. J., Witteveen, A. T., Schreiber, G. J., Kerkhoven,34825
#' R. M., Roberts, C., Linsley, P. S., Bernards, R., and Friend, S. H. (2002). Gene expression349profiling predicts clinical outcome of breast cancer.Nature, 415(6871):530-536.
#'
#' @examples
#' data(breastcancer_screened)
#' x = breastcancer_screened[,-1]
#' y = breastcancer_screened$y
#'
#' \donttest{
#' #showing distribution of correlations
#' corrs = paircorf(as.matrix(x))
#' corrplot::corrplot(corrs, method = "shade",shade.col = NA, tl.col ="black", tl.srt = 45, order = "FPC",tl.cex = 0.01)
#'
#' #showing outlier detection results via DDC
#' fit1 <- cellWise::DDC(x)
#' cellWise::cellMap(fit1$remX, fit1$stdResid, columnlabels = c(1,rep(" ",198),200), rowlabels = c(1,rep(" ",76),78), columnangle = 0,
#'         rowtitle = "Observations", columntitle = "Genes", sizetitles = 2,adjustrowlabels = 0.5, adjustcolumnlabels = 0.5)
#'
#'
#' #variable selection results via pariwise estimator and Gaussrank estimator
#' pairfit = covlasso(x,y,cor.method = "pair", scale.method = "qn", lmin = 0.1,adaptive = T)$betahat_opt
#' gaussrankfit = covlasso(x,y,cor.method = "gaussrank", scale.method = "qn", lmin = 0.1, adaptive = T)$betahat_opt
#' }
"breastcancer_screened"
