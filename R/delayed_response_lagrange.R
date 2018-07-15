#' Optimal design in delayed response situation
#'
#' @param theta Expeted effect under the alternative
#' @param sigma Standard deviation
#' @param n1 First-stage sample size
#' @param lb_n2 Lower limit for stage-two sample size
#' @param alpha Type one error rate
#' @param beta Type two error rate
#'
#' @export

lagrange_dr_design <- function(theta, sigma, n1=NA, lb_n2, alpha, beta){
  # Define inital values for integration
  N = 20
  ww <- nodes(-4, 4, N )
  h <- ( 4 - (-4) ) / ( 4 * N )
  omega = rep( 0, 4 * N + 1 )
  omega[1] = 7
  w = c( 32, 12, 32, 14 )
  omega[-1] = rep( w, N )
  omega[4*N+1] = 7



  # Compute parameters in order to obtain the same type one and two errors as JT
  th = theta / sqrt(2) / sigma
  jenturn <- jt_design(theta, sigma, 104, lb_n2, alpha, beta)

  yy <- function(x, om){
    n2 <- jenturn$n2(x)
    c2 <- jenturn$c2(x)

    y1 = om * pnorm(-c2) * dnorm(x)
    y2 = om * pnorm(-c2 + sqrt(n2) * th) * dnorm(x - sqrt(jenturn$n1) * th)

    return(c(y1, y2))
  }
  y <- apply(cbind(ww ,omega), 1, function(x) yy(x[1], x[2]))

  # Type one error
  al <- (2*h)/45 * sum(y[1,])

  # Type two error
  be <- 1 - (2*h)/45 * sum(y[2,])

  parameters <- list(alpha = al, beta = be, mu = th)


  # Define n2 function
  response_dr <- function(z, n1, cf, ce, lambda1, lambda2) {
    f <- function(v) {
      v <- v + lb_n2
      q <- 8*log( sqrt(8*pi) / ( abs(lambda2) * parameters$mu ) ) + 2 * b( parameters, z, n1, lambda1, lambda2 )
      q <- q + 4 * log( abs(v) ) + parameters$mu^2 * v + ( b( parameters, z, n1, lambda1, lambda2 ) / parameters$mu )^2 * (1/v)
      return(q)
    }
    min <- as.numeric( (-2 + sqrt( 4 + b(parameters,z,n1,lambda1,lambda2)^2 ) ) / parameters$mu^2 )
    w <- ifelse(is.finite(f(min)),f(min),1)
    if( w > 0 ){
      p <- 99999
    } else{
      ub <- min + 1
      fub <- f(ub)
      while(fub < 0){
        ub <- ub + 1
        fub <- f(ub)
      }
      p <- uniroot(f, c(min,ub))$root
    }
    p <- ifelse(z >= cf && z <= ce, p, 0)
    p <- p + lb_n2
    return(p)
  }


  # Compute the c2-function
  c2_dr <- function(z, n1, cf, ce, lambda1, lambda2){
    n2 <- round(response_dr(z, n1, cf, ce, lambda1, lambda2))
    q <- ( parameters$mu^2 * n2 - b(parameters, z, n1, lambda1, lambda2) ) / ( 2 * parameters$mu * sqrt(n2) )
    return(q)
  }



  # Compute the score
  score_dr <- function(n1, cf, ce, lambda1, lambda2){
    sc <- function(x, om){
      n <- response_dr(x, n1, cf, ce, lambda1, lambda2)
      c <- (parameters$mu^2 * n - b(parameters, x, n1, lambda1, lambda2)) / (2 * parameters$mu * sqrt(n))
      y <- (n - lb_n2) * dnorm(x - parameters$mu * sqrt(n1))
      y <- y + lambda1 * pnorm(-c) * dnorm(x)
      y <- y - lambda2 * pnorm(-c + parameters$mu * sqrt(n)) * dnorm(x - parameters$mu * sqrt(n1))
      y <- om * y
      return(y)
    }

    yy <- apply(cbind(ww, omega), 1, function(x) sc(x[1], x[2]))
    p <- (2*h)/45 * sum(yy) + n1 + lb_n2 - lambda1 * parameters$alpha + lambda2 * parameters$beta
    return(p)
  }


  # Compute lambdas and eventually n1
  fixed <- fixed(parameters)

  first_stage_dr <- function(lambda){

    if(is.na(n1)){
      optimum <- nloptr::nloptr(
        x0          = c(fixed[1]/2, 0, 2),
        eval_f      = function(x) score_dr(x[1], x[2], x[3], lambda[1], lambda[2]),
       # eval_g_ineq = function(x) err_dr(x[1], x[2], x[3], lambda[1], lambda[2]),
        lb = c(1, -3, qnorm(1-parameters$alpha)),
        ub = c(1000, qnorm(1-parameters$alpha)-0.1, 5),
        opts = list(
          algorithm = "NLOPT_LN_COBYLA",
          xtol_rel = 0.0001,
          maxeval = 999999,
          maxtime = 16200
        )
      )

      n1 <- optimum$solution[1]
      cf <- optimum$solution[2]
      ce <- optimum$solution[3]

    }else{
      optimum <- nloptr::nloptr(
        x0          = c(0, 2),
        eval_f      = function(x) score_dr(n1, x[1], x[2], lambda[1], lambda[2]),
        #eval_g_ineq = function(x) err_dr(n1, x[1], x[2], x[3], x[4]),
        lb = c(-3, qnorm(1-parameters$alpha)),
        ub = c(qnorm(1- parameters$alpha)-0.1, 5),
        opts = list(
          algorithm = "NLOPT_LN_COBYLA",
          xtol_rel = 0.0001,
          maxeval = 999999,
          maxtime = 16200
        )
      )

      cf <- optimum$solution[1]
      ce <- optimum$solution[2]
    }

    return(c(n1, cf, ce))
  }



  # Compute type one and type two error gap
  err_dr <- function(lambda1, lambda2){
    fs <- first_stage_dr(c(lambda1,lambda2))
    n1 <- fs[1]
    cf <- fs[2]
    ce <- fs[3]

    yy <- function(x, om){
      n2 <-  response_dr(x, n1, cf, ce, lambda1, lambda2)
      c2 <- (parameters$mu^2 * n2 - b(parameters, x, n1, lambda1, lambda2)) / (2 * parameters$mu * sqrt(n2))

      y1 = om * pnorm(-c2) * dnorm(x)
      y2 = om * pnorm(-c2 + sqrt(n2) * parameters$mu) * dnorm(x - sqrt(n1) * parameters$mu)

      return(c(y1, y2))
    }
    y <- apply(cbind(ww ,omega), 1, function(x) yy(x[1], x[2]))

    # Type one error
    p <- (2*h)/45 * sum(y[1,])

    # Power
    q <- (2*h)/45 * sum(y[2,])

    # Compute missing type one error resp. power
    g <- ( p - parameters$alpha )
    h <- ( 1 - parameters$beta - q )
    return(c(g, h))
  }

  if(is.na(n1)){
    l <- c(3271, 509)
  } else if(n1 == 104){
    l <- c(3125, 490)
  } else{
    l <- rootSolve::multiroot(function(x) err_dr(x[1], x[2]), c(fixed[3], fixed[4]))
  }

  # Define design

  fs <- first_stage_dr(c(l[1],l[2]))
  n1 <- fs[1]
  cf <- fs[2]
  ce <- fs[3]
  n2_out <- function(z){
    if(z >= cf && z <= ce){
      p <- response_dr(z, n1, cf, ce, l[1], l[2])
    }else{
      p <- lb_n2
    }
  }
  c2_out <- function(z){ c2_dr(z, n1, cf, ce, l[1], l[2]) }

  d <- design(cf, ce, n1, n2_out, c2_out)
  return(d)

}
