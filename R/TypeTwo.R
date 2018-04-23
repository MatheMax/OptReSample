#' Compute the type II error
#'
#' \code{type_two} computes the type II error approximately if the \code{c_2}-function and the \code{n_2}-function
#' are approximated via a set of splines.
#'
#' Note that nodes, c2 and n2 have to be of the same length.
#'
#' @param parameters Parameters (alpha, power, standardized effect) whith which you want to build your design
#' @param cf Boundary for stopping for futility after the first stage
#' @param ce Boundary for stopping for efficacy after the first stage
#' @param nodes Nodes for the stage two functions
#' @param c2 c_2-values that correspond to \code{nodes}
#' @param n1 First stage sample size
#' @param n2 n_2-values that correspond to \code{nodes}



type_two <- function(parameters, cf, ce, nodes, c2, n1, n2){
  f <- splinefun(nodes, c2)
  g <- splinefun(nodes, n2)
  N=12
  h = (ce - cf) / N
  x = seq(cf,ce,h)
  alpha=c(1,rep(2,(N-1)),1)

  y=rep(0,(N+1))
  for(i in 1:(N+1)){
    y[i] = pnorm( f(x[i]) - sqrt(abs(g(x[i])))  * parameters$mu )
    y[i] = y[i] * dnorm( x[i] - sqrt(abs(n1))  * parameters$mu )
  }
  p <- (h/2)*(t(alpha)%*%y)
  p <- p + pnorm( cf - sqrt(abs(n1)) * parameters$mu )
  return(p)
}
