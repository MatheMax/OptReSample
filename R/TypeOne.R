#' Compute the type I error
#'
#' \code{type_one} computes the type I error approximately if the \code{c_2}-function is approximated via a set of splines.
#'
#' Note that nodes and c_2 have to be of the same length.
#'
#' @param parameters Parameters (alpha, power, standardized effect) whith which you want to build your design
#' @param cf Boundary for stopping for futility after the first stage
#' @param ce Boundary for stopping for efficacy after the first stage
#' @param nodes Nodes for the c_2-function
#' @param c2 Corresponding c_2-values

type_one <- function(parameters, cf, ce, nodes, c2){
  f <- splinefun(nodes, c2)
  N = 12
  h = (ce - cf) / N
  x = seq(cf, ce,h)
  alpha=c(1, rep (2,(N-1)), 1)

  y=rep(0, (N+1))

  for(i in 1:(N+1)){
    y[i] = pnorm( f(x[i]) ) * dnorm( x[i] )
  }
  p <- (h/2)*(t(alpha)%*%y)
  p <- 1 - pnorm(cf) - p
  return(p)
}

