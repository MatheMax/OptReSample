#' Optimal design in delayed response situation
#'
#' @param theta Expeted effect under the alternative
#' @param sigma Standard deviation
#' @param n1 First-stage sample size
#' @param lb_n2 Lower limit for stage-two sample size
#' @param alpha Type one error rate
#' @param beta Type two error rate
#'
#' @return An object of class design
#'
#' @export

lagrange_dr_design <- function(theta, sigma, n1=NA, lb_n2, alpha, beta){

  # Define inital values for integration
  N = 50
  omega = rep( 0, 4 * N + 1 )
  omega[1] = 7
  w = c( 32, 12, 32, 14 )
  omega[-1] = rep( w, N )
  omega[4*N+1] = 7


  # Define parameters
  th = theta / sqrt(2) / sigma
  parameters <- list(alpha = alpha, beta = beta, mu = th)


  # Define n2 function
  response_dr <- function(z, n1, lambda1, lambda2) {
    f <- function(v) {
      v <- v + lb_n2
      q <- 4 * log(8 * pi / (lambda2 * parameters$mu)^2) + 2 * b(parameters, z, n1, lambda1, lambda2)
      q <- q + 4 * log(v) + parameters$mu^2 * v + (b(parameters, z, n1, lambda1, lambda2) / parameters$mu)^2 * (1/v)
      return(q)
    }
    min <- as.numeric( (-2 + sqrt( 4 + b(parameters,z,n1,lambda1,lambda2)^2 ) ) / parameters$mu^2 )
    w <- ifelse(is.finite(f(min)), f(min), 1)
    if( w > 0 ){
      p <- 9999999
    } else{
      ub <- min + 1
      fub <- f(ub)
      while(fub < 0){
        ub <- ub + 1
        fub <- f(ub)
      }
      p <- uniroot(f, c(min,ub))$root
    }
    p <- p + lb_n2
    return(p)
  }


  # Compute the c2-function
  c2_dr <- function(z, n1, n2, lambda1, lambda2){
    (parameters$mu^2 * n2 - b(parameters, z, n1, lambda1, lambda2)) / ( 2 * parameters$mu * sqrt(n2) )
  }


  # Compute the score
  score_dr <- function(n1, cf, ce, lambda1, lambda2){
    # Compute integral inside [cf, ce]
    ww <- nodes(cf, ce, N )
    h <- (ce - cf) / (4 * N)

    sc_int <- function(x, om){
      n <- response_dr(x, n1, lambda1, lambda2)
      c <- c2_dr(x, n1, n, lambda1, lambda2)
      y <- (n - lb_n2) * dnorm(x - parameters$mu * sqrt(n1))
      y <- y + lambda1 * pnorm(-c) * dnorm(x)
      y <- y - lambda2 * pnorm(-c + parameters$mu * sqrt(n)) * dnorm(x - parameters$mu * sqrt(n1))
      y <- om * y
      return(y)
    }

    y1 <- apply(cbind(ww, omega), 1, function(x) sc_int(x[1], x[2]))


    # Compute integrals outside this interval
    sc_out <- function(x){
      n <- lb_n2
      c <- c2_dr(x, n1, n, lambda1, lambda2)
      y <- lambda1 * pnorm(-c) * dnorm(x)
      y <- y - lambda2 * pnorm(-c + parameters$mu * sqrt(n)) * dnorm(x - parameters$mu * sqrt(n1))
      return(y)
    }

    y2 <- integrate(sc_out, -Inf, cf)$value
    y3 <- integrate(sc_out, ce, Inf)$value

    p <- (2*h)/45 * sum(y1) + y2 + y3 + n1 + lb_n2 - lambda1 * parameters$alpha + lambda2 * (1 - parameters$beta)
    return(p)
  }



  # Compute type one and type two error gap
  err_dr <- function(n1, cf, ce, lambda1, lambda2){
    # Compute integrals inside [cf, ce]
    ww <- nodes(cf, ce, N )
    h <- (ce - cf) / (4 * N)

    err_int <- function(x, om){
      n2 <- response_dr(x, n1, lambda1, lambda2)
      c2 <- c2_dr(x, n1, n2, lambda1, lambda2)

      y1 = om * pnorm(-c2) * dnorm(x)
      y2 = om * pnorm(-c2 + sqrt(n2) * parameters$mu) * dnorm(x - sqrt(n1) * parameters$mu)

      return(c(y1, y2))
    }

    y <- apply(cbind(ww, omega), 1, function(x) err_int(x[1], x[2]))


    # Compute integrals outside this interval
    to_out <- function(x){
      n2 <- lb_n2
      c2 <- c2_dr(x, n1, n2, lambda1, lambda2)
      y <- pnorm(-c2) * dnorm(x)
      return(y)
    }
    y2 <- integrate(to_out, -Inf, cf)$value + integrate(to_out, ce, Inf)$value

    tt_out <- function(x){
      n2 <- lb_n2
      c2 <- c2_dr(x, n1, n2, lambda1, lambda2)
      y <- pnorm(-c2 + sqrt(n2) * parameters$mu) * dnorm(x - sqrt(n1) * parameters$mu)
      return(y)
    }
    y3 <- integrate(tt_out, -Inf, cf)$value + integrate(tt_out, ce, Inf)$value


    # Type one error
    p <- (2*h)/45 * sum(y[1,]) + y2

    # Power
    q <- (2*h)/45 * sum(y[2,]) + y3

    # Compute missing type one error resp. power
    g <- ( p - parameters$alpha )
    h <- ( 1 - parameters$beta - q )
    return(c(g, h))
  }


  # Compute lambdas, cf, ce and eventually n1
  fixed <- fixed(parameters)

    if(is.na(n1)){
      optimum <- nloptr::nloptr(
        x0          = c(90, 0.5, 2, 3069, 483),
        eval_f      = function(x) score_dr(x[1], x[2], x[3], x[4], x[5]),
        eval_g_ineq = function(x) c( err_dr(x[1], x[2], x[3], x[4], x[5]) - c(0.001, 0.01),
                                     -err_dr(x[1], x[2], x[3], x[4], x[5]) - c(0.005, 0.01)),
        lb = c(80, 0, 1.5, 100, 10),
        ub = c(120, 1, 3, 5000, 1000),
        opts = list(
          algorithm = "NLOPT_LN_COBYLA",
          xtol_rel = 0.00001,
          maxeval = 999999,
          maxtime = 16200
        )
      )


      n1 <- round(optimum$solution[1])
      cf <- optimum$solution[2]
      ce <- optimum$solution[3]
      l1 <- ceiling(optimum$solution[4])
      l2 <- ceiling(optimum$solution[5])


    }else{
      #n1_eval <- function(n1){
      optimum <- nloptr::nloptr(
        x0          = c(0.5, 2, 3069, 483),
        eval_f      = function(x) score_dr(n1, x[1], x[2], x[3], x[4]),
        eval_g_ineq = function(x) c( err_dr(n1, x[1], x[2], x[3], x[4]) - c(0.001, 0.001),
                                     -err_dr(n1, x[1], x[2], x[3], x[4]) - c(0.001, 0)),

        lb = c(0, 1.5, 3000, 450),
        ub = c(1, 2.5, 3100, 510),
        opts = list(
          algorithm = "NLOPT_LN_COBYLA",
          xtol_rel = 0.00001,
          maxeval = 999999,
          maxtime = 16200
        )
      )

      '
      return(optimum$objective)
      }

      n1_test <- seq(80, 95, 1)
      score_test <- sapply(n1_test, n1_eval)
      opt <- n1_test[which.min(score_test)]
      '
      cf <- round(optimum$solution[1], 3)
      ce <- round(optimum$solution[2], 3)
      l1 <- ceiling(optimum$solution[3])
      l2 <- ceiling(optimum$solution[4])
    }


  while(response_dr(cf, n1, l1, l2) > 99999){
    cf <- cf + 0.001
  }

  while(response_dr(ce, n1, l1, l2) > 99999){
    ce <- ce - 0.001
  }


  n2_out <- function(z){
    if(z >= cf && z <= ce){
      p <- response_dr(z, n1, l1, l2)
    }else{
      p <- lb_n2
    }
    return(p)
  }
  c2_out <- function(z){
    n2 <- n2_out(z)
    q <- c2_dr(z, n1, n2 , l1, l2)
    return(q)
  }

  d <- design(cf, ce, n1, n2_out, c2_out)
  return(d)

}
