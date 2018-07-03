#' Evaluate the deviation of the desired error constraints
#'
#' \code{err} approximates the deviation of a design implemented with \code{lambda1} and \code{lambda2} of the error constraints that
#' were specified in \code{parameters}.
#'
#' This is a numerical approxmiation to allow a faster search for the correct lambda-values.
#' \link{find_lambda} is based on this function.
#'
#' @param parameters Parameters specifying your design
#' @param lambda1 Penalization parameter for type I error
#' @param lambda2 Penalization parameter for type II error
#'
#' @return A two-dimensional vector which gives the deviation of type I resp. type II error.

#' @export
err <- function(parameters, lambda1, lambda2) {

  n1 <- n1(parameters, lambda1, lambda2)
  c <- c_early(parameters, n1, lambda1, lambda2)
  cf <- c[1]
  ce <- c[2]


  N = 10
  h = (ce-cf)/(4*N)
  ww = nodes(cf,ce,N)

  omega = rep(0,4*N+1)
  omega[1] = 7
  wei = c(32, 12, 32, 14)
  omega[-1] = rep(wei, N)
  omega[4*N+1] = 7


  yy <- function(x, om){
    n2 <-  response(parameters, n1, lambda1, lambda2, x)
    c2 <- ( parameters$mu^2 * n2 - b(parameters, x, n1, lambda1, lambda2) ) / ( 2 * parameters$mu * sqrt(abs(n2)) )

    y1 = om * pnorm( c2 ) * dnorm( x )
    y2 = om * pnorm( c2 - sqrt(abs(n2)) * parameters$mu) * dnorm( x - sqrt(n1) * parameters$mu )

    return(c(y1,y2))
  }

  y <- sapply(ww, yy)

  q <- (2*h)/45 * sum(y[2,])
  q <- 1 - pnorm ( cf - sqrt( n1 ) * parameters$mu ) - q

  p <- (2*h)/45 * sum(y[1,])
  p <- 1 - pnorm(cf) - p

  g <- ( p - parameters$alpha )
  h <- ( 1 - parameters$beta - q )
  return(c(g,h))
}


#' Starting values for lambda
#'
#' \code{lambda_start} gives suitable starting values to find the values of lambda that yield to desired
#' type I - and type II - error.
#'
#' These values are used as starting values in \link{find_lambda}.
#' The choice of the starting values seems quite arbitrarily. It is based on own examples. Improvments are very welcome.
#'
#' @param parameters Parameters specifying the design.
#' @export


lambda_start <- function(parameters){
  if( parameters$mu>0.5 ){
    l2 <- 25
  } else if( parameters$mu>0.4 ){
    l2 <- 2000 * parameters$alpha
  } else if( parameters$mu>0.3 ){
    l2 <- max( 150 - 1500 * parameters$alpha, 1)
  } else if( parameters$mu>0.25 ){
    l2 <- max( 200 - 1500 * parameters$alpha, 1)
  } else{
    l2 <- max( 300 * ( 1 + (0.2 - parameters$mu) * 10 )^2, 1)
  }

  ratio <- parameters$beta / parameters$alpha

  if(ratio<=9){
    g <- function(z) {
      x = c(1,2,3,4,6,8,9)
      y = c(2,2.5,3,3.4,4,6,6.5)
      q <- lm(y~x)
      p <- q$coefficients[1] + q$coefficients[2] * z
      return(p)
    }
    l1 <- as.numeric( g(ratio) * l2 )
  }  else{
    g <- function(z) {
      x = c(10, 11, 12, 13)
      y = c(11, 10.5, 6.6, 6.6)
      q <- lm(y ~ x + I(x^2))
      p <- q$coefficients[1] + q$coefficients[2] * z + q$coefficients[3] * z^2
      return(p)
    }
    l1 <- as.numeric( g(ratio) * l2 )
  }
  l1 <- max(l1, 2 * l2)

  return(c(l1, l2))

}


#' Find the correct values of lambda
#'
#' \code{find_lambda} returns the values of lambda that yield to the desired type I - and type II - error.
#'
#' Use these values in \link{lagrange_design} to obtain an optimal design with type I - and type II - error rate
#' as prespecified.
#'
#' @param parameters Parameters specifying the design.
#'
#' @return A two-dimensional vector (lambda1, lambda2).
#' @export

find_lambda <- function(parameters){
  l <- lambda_start(parameters)

  err_opt<-function(lambda1, lambda2){err(parameters, lambda1, lambda2)}
  opt <- rootSolve::multiroot(function(x) err_opt(x[1], x[2]), c(l[1], l[2]), atol=c(0.01,0.01))

  return(ceiling(opt$root))
}


#' Find lambda values faster
#'
#' \code{find_lambda_direct} is a way to find approximately correct lambda-values faster than with \link{find_lambda}.
#'
#' This function is based on \link{direct_design}. Note therefore that is does not yield a purely Lagrangian approach
#' anymore and is more liable to numerical errors.
#'
#' @param parameters Parameters specifying the design.
#'
#' @return A two-dimensional vector (lambda1, lambda2).
#' @export


find_lambda_direct <- function(parameters){
  d <- direct_design(parameters)

  l <- lambda_start(parameters)

  f <- function(lambda1,lambda2){
    c <- c_early(parameters, d$n1, lambda1, lambda2)
    p = rep(0,2)
    p[1] = d$cf - c[1]
    p[2] = d$ce - c[2]
    return(p)
    }

  opt <- rootSolve::multiroot(function(x) f(x[1], x[2]), c(l[1], l[2]))

  return(opt$root)

}


