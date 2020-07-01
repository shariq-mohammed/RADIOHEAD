#' FDR-based variable selection
#'
#' Computes the estimated coefficients of \eqn{\beta} from the radiogenomics model
#'
#' @param b MCMC samples of the parameter \eqn{\beta} - dimension \eqn{N} samples x \eqn{L} predictors. Output from the groupSS function
#' @param x_cnames names of columns in the predictor matrix
#' @param ... other arguments to \eqn{localFDR()}
#' @keywords fdr_var_selection()
#' @export
#'
#' @importFrom stats density
#' @return est_coef: estimated coefficients after the FDR-based variable selection
#' @examples
#' b = matrix(rnorm(100*5), nrow=100)
#' x_cnames = paste0('Col',1:5)
#' fdr_var_selection(b,x_cnames)

####################################
### FDR-based variable selection ###
####################################

fdr_var_selection = function(b, x_cnames, ... ){
  post_inc_prob = 1-localFDR(b,...)

  if(sum(!is.na(post_inc_prob))>0){
    inds = which(!is.na(post_inc_prob))
    est_coef = apply(b, 2,
                     function(x){
                       obj = density(x)
                       m = obj$x[which.max(obj$y)]
                       if(m == max(obj$x[obj$x<0]) | m == min(obj$x[obj$x>0])) return(0)
                       else return(m)
                     })[inds]
    names(est_coef) = x_cnames[inds]
    return(est_coef)
  }else return('No significant associations.')
}
