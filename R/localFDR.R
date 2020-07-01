#' Local False Discovery Rates
#'
#' Computes the local false discovery rates (FDR) for the paramters given MCMC samples for those parameters
#'
#' @param chain MCMC samples of the parameter - dimension \eqn{N} samples x \eqn{L} parameters.
#' @param thres threshold used to compute the Local FDR. Defaults to 0.001.
#' @param FDRc level at which the FDR is controlled (Bayesian FDR). Defaults to 0.05.
#' @keywords localFDR()
#' @export
#' @return r: local FDRs of the \eqn{L} predictors
#' @examples
#' chain = matrix(rnorm(100*5), nrow=100)
#' localFDR(chain)

###################################
### Local False Discovery Rates ###
###################################
localFDR = function(chain, thres = 0.001, FDRc = 0.05){
  mIter = nrow(chain)
  p = ncol(chain)
  x = sapply(1:p, function(k) mean(abs(chain[,k])<thres))
  x2 = cumsum(x[order(x)])/1:p
  r = rep(NA, p)
  x.ind = order(x)[which(as.numeric(x2<FDRc)==1)]
  r[x.ind] = x[x.ind]
  r
}
