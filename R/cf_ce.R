#' d Score / d c_f
#'
#' \code{d_cf} gives the gradient of the score w.r.t. c_f.  Mainly needed for \link{c_early}
#' @param parameters The parameters (alpha, power, etc.) you want to use for your design
#' @param n1 The first-stage sample size
#' @param lambda1 The penalization parameter for type one error
#' @param lambda2 The penalization parameter for type two error
#' @param x The z-value of the first stage
d_cf <- function(parameters,n1,lambda1,lambda2,x){
  #n2 berechnen
  n2 <- response(parameters,n1,lambda1,lambda2,x)

  #c2 berechnen
  c2 <- ( a(parameters)^2 * n2 - b(parameters,x,n1,lambda1,lambda2) ) / ( 2 * a(parameters) * sqrt(abs(n2)) )

  p <- - lambda1 * dnorm(x) + lambda2 * dnorm( x - a(parameters) * sqrt(abs(n1)) )
  p <- p - n2 * dnorm( x - sqrt(abs(n1)) * a(parameters) )
  p <- p + lambda1 * pnorm(c2) * dnorm(x)
  p <- p - lambda2 * pnorm(c2 - a(parameters) * sqrt(abs(n2)) ) * dnorm( x - a(parameters) * sqrt(abs(n1)) )
  return(p)
}


#' d Score / d c_e
#'
#' \code{d_ce} gives the gradient of the score w.r.t. c_f.  Mainly needed for \link{c_early}
#' @param parameters The parameters (alpha, power, etc.) you want to use for your design
#' @param n1 The first-stage sample size
#' @param lambda1 The penalization parameter for type one error
#' @param lambda2 The penalization parameter for type two error
#' @param x The z-value of the first stage
d_ce <- function(parameters,n1,lambda1,lambda2,x){
  #n2 berechnen
  n2 <- response(parameters,n1,lambda1,lambda2,x)

  #c2 berechnen
  c2 <- ( a(parameters)^2 * n2 - b(parameters,x,n1,lambda1,lambda2) ) / ( 2 * a(parameters) * sqrt(abs(n2)) )

  p <- n2 * dnorm( x - sqrt(abs(n1)) * a(parameters) )
  p <- p - lambda1 * pnorm(c2) * dnorm(x)
  p <- p + lambda2 * pnorm(c2 - a(parameters) * sqrt(abs(n2)) ) * dnorm( x - a(parameters) * sqrt(abs(n1)) )
  return(p)
}


#' Optimal Stopping Boundaries
#'
#' \code{c_early} gives the optimal boundaries for early stopping of an adaptive design
#'
#' The parameter \code{n1} has to be given as well as the penalty parameters \code{lambda_1} and \code{lambda_2}.
#' The boundaries are computed as roots of the gradient of the score w.r.t. c_f and c_e.
#'
#' @param parameters The parameters (alpha, power, etc.) you want to use for your design
#' @param n1 The first-stage sample size
#' @param lambda1 The penalization parameter for type one error
#' @param lambda2 The penalization parameter for type two error
#' @return A two dimensional vector (c_f, c_e)


c_early <- function(parameters,n1,lambda1,lambda2){
  lb <- qnorm(1-parameters$alpha)
  z <- d_ce(parameters,n1,lambda1,lambda2,lb)
  ub <- lb + 0.1
  z2 <- d_ce(parameters,n1,lambda1,lambda2,ub)

  while(sign(z)==sign(z2)){
    ub <- ub + 0.1
    z2 <- d_ce(parameters,n1,lambda1,lambda2,ub)
  }

  ce <- uniroot(function(x){d_ce(parameters,n1,lambda1,lambda2,x)},c(lb,ub),extendInt="no")$root


  lb2 <- ce - 0.1
  z3 <-  d_cf(parameters,n1,lambda1,lambda2,ce)
  z4 <-  d_cf(parameters,n1,lambda1,lambda2,lb2)

  while(sign(z3) == sign(z4)){
    lb2 <- lb2 - 0.1
    z4 <- d_cf(parameters,n1,lambda1,lambda2,lb2)
  }

  cf <- uniroot(function(x){d_cf(parameters,n1,lambda1,lambda2,x)},c(lb2,ce),extendInt="no")$root


  return(c(cf,ce))
}
