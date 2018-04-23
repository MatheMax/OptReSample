#' Standardized effect size
#'
#' For non-standardized effects, \code{a} returns the standardized effect size
#'
#' @param parameters Parameters specifying the design

a <- function(parameters) {
   as.numeric( ( parameters$mualt - parameters$mu0 ) / parameters$sigma )
}

#'The value b as defined in the manuscript
#'
#'\code{b} is a help function to compute optimal designs in the Lagrangian framework
#'
#' @param parameters Parameters specifying the design
#' @param z First stage z value
#' @param n1 First stage sample size
#' @param lambda1 Penalization parameter for type I error
#' @param lambda2 Penalization parameter for type II error

b <- function(parameters,z,n1,lambda1,lambda2) {
   2 *  log(abs(lambda2 / lambda1)) + 2 * sqrt(n1) * a(parameters) * z - n1 * a(parameters)^2
}


#' Define equi distance nodes
#'
#' \code{nodes} defines equi distance nodes inside the invterval (l,u)
#'
#' @param l Lower Bound
#' @param u Uper Bound
#' @param N In total, there are 4N+1 nodes

nodes <- function(l,u,N) {
  h=(u-l)/(4*N)
  x=seq(l,u,h)
  return(x)
}

