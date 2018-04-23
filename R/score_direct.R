#' Compute the score
#'
#' \code{score_direct} computes the objective score approximately if the \code{c_2}-function and the \code{n_2}-function
#' are approximated via a set of splines.
#'
#' Note that nodes, c2 and n2 have to be of the same length.
#' The score is given as the expected sample size of the design under the alternative hypothesis.
#' The integrals are approximated via the trapezial rule.
#'
#'
#' @param parameters Parameters (alpha, power, standardized effect) whith which you want to build your design
#' @param cf Boundary for stopping for futility after the first stage
#' @param ce Boundary for stopping for efficacy after the first stage
#' @param nodes Nodes for the stage two functions
#' @param c2 c_2-values that correspond to \code{nodes}
#' @param n1 First stage sample size
#' @param n2 n_2-values that correspond to \code{nodes}




score_direct <- function(parameters,cf,ce,nodes,c2,n1,n2){
  g <- splinefun(nodes,n2)

  N=12

  h = (ce - cf) / N
  x = seq(cf,ce,h)
  alpha=c(1,rep(2,(N-1)),1)

  y=rep(0,(N+1))
  for(i in 1:(N+1)){
    alpha[i] = alpha[i] * g(x[i])
    y[i] <- dnorm( x[i] - sqrt(abs(n1)) *a(parameters) )
  }
  p <- (h/2)*(t(alpha)%*%y)
  p <- p + n1
  return(p)
}

