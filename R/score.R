#' Score in the Lagrangian framework
#'
#' \code{score} defines the score in the Lagrangian framework.
#'
#' The score is given as the expected sample size of the trial under the alternative hypothesis.
#' @param parameters Parameters specifying the design
#' @param cf Boundary for stopping for futility
#' @param ce Boundary for stopping for efficacy
#' @param n1 First stage sample size
#' @param lambda1 Penalization parameter for type I error
#' @param lambda2 Penalization parameter for type II error
#'
#' @export

score <- function( parameters, cf, ce, n1, lambda1, lambda2 ) {
  # Compute the integral
  N = 20
  x <- nodes( cf, ce, N )
  h <- ( ce - cf ) / ( 4 * N )
  omega = rep( 0, 4 * N + 1 )
  omega[1] = 7
  w = c( 32, 12, 32, 14 )
  omega[-1] = rep( w, N )
  omega[4*N+1] = 7

  sc <- function(x) {
    n <- response( parameters, n1, lambda1, lambda2, x )
    c <- ( parameters$mu^2 * n - b( parameters, x, n1, lambda1, lambda2) ) / (2 * parameters$mu * sqrt(n))
    y <- n * dnorm( x - parameters$mu * sqrt(n1) )
    y <- y - ( lambda1 * pnorm(c) * dnorm(x) )
    y <- y + ( lambda2 * pnorm( c - parameters$mu * sqrt(n) ) * dnorm( x - parameters$mu * sqrt(n1)))
    return(y)
  }


  y <- sapply(x, sc)

  p <- (2 * h) / 45 * (t(omega) %*% y)


  # Compute the score
  g <- n1 + lambda1 * ( 1 - parameters$alpha - pnorm(cf) )
  g <- g + lambda2 * ( pnorm( cf - parameters$mu * sqrt(n1) ) - parameters$beta )
  g <- g + p

  return(g)
}
