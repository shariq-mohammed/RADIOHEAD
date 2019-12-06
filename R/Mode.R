#' Mode for a vector of numbers
#'
#' Computes mode for a given set of MCMC samples of a parameter
#'
#' @param x numeric vector of \eqn{N} MCMC samples (or any numeric vector)
#' @keywords Mode()
#' @export
#' @return m mode for the group of numbers
#' @examples Mode()

####################################
### Mode for a vector of numbers ###
####################################

Mode = function(x){
  obj = density(x)
  m = obj$x[which.max(obj$y)]
  if(m == max(obj$x[obj$x<0]) | m == min(obj$x[obj$x>0])) return(0)
  else return(m)
}