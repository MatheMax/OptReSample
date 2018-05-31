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

  sc <- function(x, parameters_=parameters, cf_=cf, ce_=ce, n1_=n1, lambda1_=lambda1, lambda2_=lambda2) {
    n <- response( parameters_, n1_, lambda1_, lambda2_, x )
    c <- ( parameters_$mu^2 * n - b( parameters_, x, n1_, lambda1_, lambda2_ ) ) / ( 2 * parameters_$mu * sqrt(n) )
    y <- n * dnorm( x - parameters_$mu * sqrt(n1_) )
    y <- y - ( lambda1_ * pnorm(c) * dnorm(x) )
    y <- y + ( lambda2_ * pnorm( c - parameters_$mu * sqrt(n) ) * dnorm( x - parameters_$mu * sqrt(n1_) ) )
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
