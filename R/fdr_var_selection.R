#' FDR-based variable selection
#'
#' Computes the estimated coefficients of \eqn{\beta} from the radiogenomics model
#'
#' @param b MCMC samples of the parameter \eqn{\beta} - dimension \eqn{N} samples x \eqn{L} predictors. Output from the groupSS function
#' @param x_cnames names of columns in the predictor matrix
#' @param ... other arguments to \eqn{pep()}
#' @keywords fdr_var_selection()
#' @export
#' @return est_coef estimated coefficients after the FDR-based variable selection 
#' @examples fdr_var_selection()

####################################
### FDR-based variable selection ###
####################################

fdr_var_selection = function(b,
                             x_cnames,
                             ...
                             ){
  post_inc_prob = 1-RadioGenLGG:::pep(b,...)
  
  if(sum(!is.na(post_inc_prob))>0){
    inds = which(!is.na(post_inc_prob))
    est_coef = apply(b, 2, RadioGenLGG:::Mode)[inds]
    names(est_coef) = x_cnames[inds]
    return(est_coef)
  }else return('No significant associations.')
}