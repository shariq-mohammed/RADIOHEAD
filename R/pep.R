#' Posterior Exclusion Probabilities
#'
#' Computes the posterior exclusion probabilities for the paramters given MCMC samples for those parameters
#'
#' @param chain MCMC samples of the parameter - dimension \eqn{N} samples x \eqn{L} parameters.
#' @param thres threshold used to compute the posterior exclusion probabilities. Defaults to 0.001.
#' @param FDRc level at which the FDR is controlled (average Bayesian FDR). Defaults to 0.05.
#' @keywords pep()
#' @return r posterior exclusion probabilities of the \eqn{L} predictors
#' @examples pep()

#########################################
### Posterior Exclusion Probabilities ###
#########################################
pep = function(chain,
               thres = 0.001,
               FDRc = 0.05){
  mIter = nrow(chain)
  p = ncol(chain)
  x = sapply(1:p, function(k) mean(abs(chain[,k])<thres))
  x2 = cumsum(x[order(x)])/1:p
  r = rep(NA, p)
  x.ind = order(x)[which(as.numeric(x2<FDRc)==1)]
  r[x.ind] = x[x.ind]
  r
}
