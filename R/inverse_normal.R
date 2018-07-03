#' Compute optimal sample size based on inverse normal method
#'
#' @export

inverse_normal_design <- function(parameters){
  n_fixed <- fixed(parameters)[1]
  cf <- 0
  ce <- 4
  h = (ce - cf) / 10
  ww = seq(cf, ce, h)
  alph = c(1, rep(2, 9), 1)


  c2_in <- function(z, n1){
    w <- rep(0,2)
    w[1] <- sqrt(n1 / n_fixed)
    w[2] <- sqrt( 1 - w[1]^2 )
    p <- (qnorm(1 - parameters$alpha) - w[1]*z)/w[2]
    return(p)
  }

  n2_in <- function(z, n1, lambda){
    c <- c2_in(z, n1)
    f <- function(n){
      log(8 * pi / (lambda * parameters$mu)^2) + log(n) + (c - sqrt(n) * parameters$mu)^2
    }
    w <- f(1)
    if( w > 0 ){
      p <- 0
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

  '
  # Type one error
  to_in <- function(n1, cf, ce){
    p <- 1 - pnorm(cf)
    h = (ce - cf) / 10
    ww = seq(cf, ce, h)
    alph = c(1, rep(2, 9), 1)

    # z <- (z, alpha)
    f <- function(z, al){ al * pnorm(c2_in(z, n1)) * dnorm(z) }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <-  p - (h/2) * sum(y)
    return(p)

  }

  power_lack <- function(n1, cf, ce, lambda){
    h = (ce - cf) / 10
    ww = seq(cf, ce, h)

    # Define type two error
    f <- function(z, al){
      al * pnorm(c2_in(z, n1) - sqrt(n2_in(z, n1, lambda)) * parameters$mu) * dnorm(z - sqrt(n1) * parameters$mu)
    }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <- parameters$beta - (h/2) * sum(y) - pnorm(cf - sqrt(n1) * parameters$mu)
    return(p)
  }
  '


  score_in <- function(n1, lambda){
    # z <- (z, alpha)
    f <- function(z, al){
      n2 <- round(n2_in(z, n1, lambda))
            q <- al * n2 * dnorm(z - sqrt(n1) * parameters$mu) +
           al * lambda * pnorm(c2_in(z, n1) - sqrt(n2) * parameters$mu) *
                    dnorm(z - sqrt(n1) * parameters$mu)
      return(q)
    }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <- (h/2) * sum(y)
    p <- p + n1 + lambda * ( pnorm(cf - sqrt(n1) * parameters$mu) - parameters$beta )
    return(p)
  }

  n1_in <- function(lambda){
    optimum <- nloptr::nloptr(
      x0          = c(n_fixed/2),
      eval_f      = function(x) score_in(x, lambda),
      lb = c(0.1 * n_fixed),
      ub = c(0.9 * n_fixed),
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
      n1 <- n1_in(lambda)

      # Define type two error
      f <- function(z, al){
        al * pnorm(c2_in(z, n1) - sqrt(n2_in(z, n1, lambda)) * parameters$mu) * dnorm(z - sqrt(n1) * parameters$mu)
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

  n1_out <- round(n1_in(lambda_opt))

  n2_out2 <- function(z){ n2_in(z, n1_out, lambda_opt) }
  c2_out <- function(z){ c2_in(z, n1_out) }

  xx <- seq(cf, ce, 0.01)
  n2_test <- sapply(xx, n2_out2)
  cf_out <-  max( xx[min(which(n2_test > 3))] , cf)
  ce_out <-  min( xx[max(which(n2_test > 3))], ce)

  n2_out <- function(z){
    ifelse(cf_out <= z && z <= ce_out, n2_out2(z), 0)
  }
  d <- design(cf_out, ce_out, n1_out, n2_out, c2_out)
  return(d)

}


