#'A simplified version of the function b that is defined in the manuscript
#'
#'\code{b} is a help function to compute optimal designs in the Lagrangian framework
#'
#' @param parameters Parameters specifying the design
#' @param z First stage z value
#' @param n1 First stage sample size
#' @param lambda1 Penalization parameter for type I error
#' @param lambda2 Penalization parameter for type II error
#' @export

b <- function(parameters, z, n1, lambda1, lambda2) {
   2 *  log( abs(lambda2 / lambda1) ) + 2 * sqrt(n1) * parameters$mu * z - n1 * parameters$mu^2
}


#' Define equi distance nodes
#'
#' \code{nodes} defines equi distance nodes inside the invterval (l,u)
#'
#' @param l Lower Bound
#' @param u Uper Bound
#' @param N In total, there are 4N+1 nodes
#' @export

nodes <- function(l, u, N) {
  h=(u-l)/(4*N)
  x=seq(l,u,h)
  return(x)
}


#' Define weights for Milne integration
#'
#' @param N 4*N+1 nodes
#'
#' @export
#'

omega <- function(N){
  omega = rep(0, 4*N+1)
  omega[1] = 7
  se = c(32, 12, 32, 14)
  omega[-1] = rep(se, N)
  omega[4*N+1] = 7
  return(omega)
}



