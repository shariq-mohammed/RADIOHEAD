#' Geometric Principal Component Analysis
#'
#' Computes the principal component (PC) scores from the geometric principal
#' component analysis on probability density functions.
#'
#' @param x a list with each entry as a vector of values required to
#'          compute the probability density function for that case/subject
#' @param perc.var percentage variance explained by PCs (between 0 to 100)
#' @param m number of equally spaced points at which the densities are to be estimated
#'
#' @return pcScores principal component scores; a matrix of size \eqn{n}x\eqn{p},
#'                  where \eqn{p} is the number of PCs included
#'
#' @export
#'
#' @examples
#' x = lapply(1:10, function(i) rnorm(100))
#' PCScores = geomPCA(x)

geomPCA = function(x, perc.var = 99, m = 512){
  perc.var = perc.var/100
  stopifnot(perc.var > 0)
  if(perc.var>1) perc.var = 1

  n = length(x)

  min_coef = min(unlist(x), na.rm=T)
  max_coef = max(unlist(x), na.rm=T)
  den = matrix(nrow = n, ncol = m)
  for(i in 1:n){
    temp_vec = c(x[[i]][!is.na(x[[i]])])
    den[i,] = density(temp_vec, from = min_coef, to = max_coef,
                      n = m, bw='nrd')$y
    den[i,] = den[i,]/(max_coef - min_coef)
    den[i,] = den[i,]/sum(diff(seq(0,1,length = m))*(den[i,-m]+ den[i,-1])/2)
  }

  sqrt_den = sqrt(den)
  eps = .5
  nmv = 100
  kmean_sqrt = rep(1, m)

  iter = 1
  while(nmv[iter] > 1e-10){
    vbar = rep(0, m)
    v = matrix(nrow = n, ncol = m)
    for(i in 1:n){
      if(abs(sum(kmean_sqrt-(sqrt_den[i,])))<1e-10){
        v[i,] = rep(0,m)
      } else{
        th = acos(sum(diff(seq(0,1,length=m))*
                        ((kmean_sqrt*(sqrt_den[i,]))[-m]+
                           (kmean_sqrt*(sqrt_den[i,]))[-1])/2))
        v[i,] = (th/sin(th))*((sqrt_den[i,])-kmean_sqrt*cos(th))
      }

      vbar = vbar + v[i,]/n
    }

    nv = sqrt(sum((((eps*vbar)*(eps*vbar))[-m]+
                     ((eps*vbar)*(eps*vbar))[-1])/2)/(m-1))
    kmean_sqrt = cos(nv)*kmean_sqrt + sin(nv)*(eps*vbar)/nv
    iter = iter+1
    nmv = append(nmv, sqrt(sum(((vbar*vbar)[-m]+(vbar*vbar)[-1])/2)/(m-1)))
  }

  ieden = matrix(0, nrow = n, ncol = m)
  for(i in 1:n){
    if(abs(sum(kmean_sqrt-(sqrt_den[i,])))<1e-10){
      ieden[i,] = rep(0,m)
    } else{
      th = acos(sum(diff(seq(0,1,length=m))*
                      ((kmean_sqrt*(sqrt_den[i,]))[-m]+
                         (kmean_sqrt*(sqrt_den[i,]))[-1])/2))
      ieden[i,] = (th/sin(th))*((sqrt_den[i,])-kmean_sqrt*cos(th))
    }
  }

  K = matrix(0, nrow = m, ncol = m)
  for(i in 1:n) K = K+tcrossprod(ieden[i,])
  K = K/(n-1)

  svdecomp = svd(K)
  p = sum(cumsum((svdecomp$d^2)/sum(svdecomp$d^2))<perc.var)+1
  U = svdecomp$u[,1:p]
  pcScores = ieden%*%U

  pcScores
}
