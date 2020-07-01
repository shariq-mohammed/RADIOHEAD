#' MCMC algorithm for the estimation of the Bayesian model
#'
#' Gibbs sampling algorithm for the estiamtion of the model with group spike-and-slab prior
#'
#' @param y response for \eqn{n} subjects
#' @param X predictors - dimension \eqn{n} subjects x \eqn{L} predictors
#' @param g_id numeric vector indicating group membership in the predictors. Contains values \eqn{1,\ldots,g} if the predictors belong to \eqn{g} groups
#' @param beta_init initial values for \eqn{\beta}. Defaults to 0.
#' @param zeta_init initial values for \eqn{\zeta}. Defaults to random generation from Ber(0.5).
#' @param w_init initial value for \eqn{w} - the complexity parameter. Defaults to proportion of 1s in zeta_init.
#' @param nusq_init initial values for \eqn{\nu^2}. Defaults to 1.
#' @param sigsq_init initial values for \eqn{\sig^2}. Defaults to 1.
#' @param v0 hyperparameter \eqn{v0}. Defaults to 0.005.
#' @param a1 hyperparameter \eqn{a1}. Defaults to 0.001.
#' @param a2 hyperparameter \eqn{a2}. Defaults to 0.001.
#' @param c1 hyperparameter \eqn{c1}. Defaults to 0.001.
#' @param c2 hyperparameter \eqn{c2}. Defaults to 0.001.
#' @param Nmcmc number of MCMC samples to generate. Defaults to 5000.
#' @param ind indices of MCMC samples to use after burnin and thinning. Defaults to 1 to Nmcmc.
#' @keywords groupSS()
#' @export
#'
#' @importFrom stats rbeta rbinom rgamma rnorm
#' @return zeta: MCMC samples of \eqn{\zeta} - dimension \eqn{N} samples x \eqn{g} groups
#' @return b: MCMC samples of \eqn{\beta} - dimension \eqn{N} samples x \eqn{L} predictors
#' @return x_cnames: column names from the predictor matrix
#' @examples
#' p = 10
#' n = 100
#' beta = rep(c(2,0), c(p/2,p/2))
#' g_id = rep(c(1,2), c(p/2,p/2))
#' X = matrix(rnorm(n*p), nrow = n)
#' y = rnorm(n, mean = X%*%beta, sd = 0.01)
#' groupSS(y, X, g_id, Nmcmc = 50)

###############################################################
### MCMC algorithm for the estimation of the Bayesian model ###
###############################################################

groupSS = function(y, X, g_id, beta_init = NULL, zeta_init = NULL, w_init = NULL,
                   nusq_init = NULL, sigsq_init = 1, v0 = 0.005, a1 = 0.001,
                   a2 = 0.001, c1 = 0.001, c2 = 0.001, Nmcmc = 5000, ind = 1:Nmcmc){
  nc = ncol(X)
  nr = nrow(X)
  stopifnot(is.numeric(g_id), min(g_id)>=1, g_id%%1==0)
  g = max(g_id)

  if(is.null(beta_init)) beta_init = t(rep(0,nc))
  if(is.null(zeta_init)){
    zeta_init = rbinom(g,1,0.5)
    zeta_init[zeta_init==0] = v0
  } else stopifnot(length(zeta_init)==g, unique(zeta_init) %in% c(1,v0))
  if(is.null(w_init)) w_init = sum(zeta_init[g_id]==1)/nc
  if(is.null(nusq_init)) nusq_init = rep(1, nc)

  b = beta_init
  zeta = zeta_init
  w = w_init
  nusq = nusq_init
  sigsq = sigsq_init

  zeta_store = matrix(nrow = Nmcmc, ncol = g)
  b_store = matrix(nrow = Nmcmc, ncol = nc)

  xtx = crossprod(X)
  xty = crossprod(X, y)
  GMinv = diag(1/(nusq*zeta[g_id]))

  for(i in 1:Nmcmc){
    if(i %% 1000 == 0) print(paste(Sys.time(),'Iteration',i,sep=" "))

    bVar.chol = chol(xtx+GMinv)
    bVar = chol2inv(bVar.chol)
    bMean = crossprod(bVar, xty)

    b = as.vector(backsolve(bVar.chol/sqrt(sigsq), rnorm(nc))) + bMean
    b_store[i,] = b

    tempExp = exp(-tapply((b^2)/(2*sigsq*nusq),g_id,sum))
    pr1 = (1-w)*(tempExp^(1/v0))/(sqrt(v0)^table(g_id))
    pr2 = w*tempExp
    temp.probs = pr2/(pr1+pr2)
    temp.probs[is.na(temp.probs)] = 0.5 #case when pr1 and pr2 are extremely small
    zeta = rbinom(g, 1, prob = temp.probs)
    zeta[zeta==0] = v0
    zeta_store[i,] = zeta

    nusq = 1/rgamma(nc, shape = a1+0.5, rate = a2+((b^2)/(2*sigsq*zeta[g_id])))

    w = rbeta(1, 1+sum(zeta==1), 1+sum(zeta==v0))

    GMinv = diag(1/(nusq*zeta[g_id]))
    sigsq = 1/rgamma(1, shape = c1+((nr+nc)/2),
                     rate = c2+((sum((y-crossprod(t(X),b))^2)+
                                   sum((b^2)*diag(GMinv)))/2))
  }

  list(zeta = zeta_store[ind,],
       b = b_store[ind,],
       x_cnames = colnames(X)
  )
}
