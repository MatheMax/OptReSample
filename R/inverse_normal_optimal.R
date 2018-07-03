#' Compute the optimal design based on the inverse normal combination test
#'
#' @export

optimal_inverse_normal_design <- function(parameters){
  n_fixed <- fixed(parameters)[1]

  c2_in <- function(z, w){
    c <- qnorm(1 - parameters$alpha)
    p <- (c - w * z) / sqrt( 1 - w^2 )
    return(p)
  }

  n2_in <- function(z, n1, cf, ce, w, lambda){
    c <- c2_in(z, w)
    f <- function(n){
      log(8 * pi / (lambda * parameters$mu)^2) + log(n) + (c - sqrt(n) * parameters$mu)^2
    }
    w <- f(1)
    if( w > 0 ){
      p <- 99999
    } else{
      ub <- 5
      fub <- f(ub)
      while(fub<0){
        ub <- ub + 5
        fub <- f(ub)
      }
      p <- uniroot(f,c(1,ub))$root
    }

    return(p)
  }



  score_in <- function(n1, cf, ce, w, lambda){
    h = (ce - cf) / 10
    ww = seq(cf, ce, h)
    alph = c(1, rep(2, 9), 1)

    # z <- (z, alpha)
    f <- function(z, al){
      n2 <- round(n2_in(z, n1, cf, ce, w, lambda))
      q <- al * n2 * dnorm(z - sqrt(n1) * parameters$mu) +
        al * lambda * pnorm(c2_in(z, w) - sqrt(n2) * parameters$mu) *
        dnorm(z - sqrt(n1) * parameters$mu)
      return(q)
    }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <- (h/2) * sum(y)
    p <- p + n1 + lambda * ( pnorm(cf - sqrt(n1) * parameters$mu) - parameters$beta )
    return(p)
  }

  to_in <- function(n1, cf, ce, w){
    p <- 1 - pnorm(cf)
    h = (ce - cf) / 10
    ww = seq(cf, ce, h)
    alph = c(1, rep(2, 9), 1)

    # z <- (z, alpha)
    f <- function(z, al){ al * pnorm(c2_in(z, w)) * dnorm(z) }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <-  p - (h/2) * sum(y)
    return(p)

  }


  first_stage_in <- function(lambda){
    optimum <- nloptr::nloptr(
      x0          = c(n_fixed/2, 0, 2, 0.5),
      eval_f      = function(x) score_in(x[1], x[2], x[3], x[4], lambda),
      eval_g_ineq = function(x) to_in(x[1], x[2], x[3], x[4]) - parameters$alpha,
      lb = c(0.1 * n_fixed, -1, qnorm(1-parameters$alpha), 0.001),
      ub = c(0.9 * n_fixed, qnorm(1 - parameters$alpha)-0.1, 4, 0.99),
      opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 0.0001,
        maxeval = 999999,
        maxtime = 16200
      )
    )

    return(optimum$solution)
  }

  # Find lambda s.t. power constraint holds
  optimal_lambda <- function(parameters){

    power_lack <- function(lambda){
      fs <- first_stage_in(lambda)
      n1 <- fs[1]
      cf <- fs[2]
      ce <- fs[3]
      w <- fs[4]

      h = (ce - cf) / 10
      ww = seq(cf, ce, h)
      alph = c(1, rep(2, 9), 1)



      # Define type two error
      f <- function(z, al){
        al * pnorm(c2_in(z, w) - sqrt(n2_in(z, n1, cf, ce, w, lambda)) * parameters$mu) *
          dnorm(z - sqrt(n1) * parameters$mu)
      }
      y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
      p <- parameters$beta - (h/2) * sum(y) - pnorm(cf - sqrt(n1) * parameters$mu)
      return(p)
    }
    lam <- try(uniroot(power_lack, c(10, 20000), extendInt="yes")$root, silent=T)
    if(class(lam)=="try-error"){ lam = 999999 }
    return(as.numeric(lam))
  }

  lambda_opt <- ceiling(optimal_lambda(parameters))

  fs <- first_stage_in(lambda_opt)
  n1 <- round(fs[1])
  cf <- fs[2]
  ce <- fs[3]
  w <- fs[4]


  n2_out <- function(z){ n2_in(z, n1, cf, ce, w, lambda_opt) }
  c2_out <- function(z){ c2_in(z, w) }

  d <- design(cf, ce, n1, n2_out, c2_out)
  return(d)

}
