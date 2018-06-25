#' Defines design based on conditional power
#'
#' @export

cond_power_design <- function(parameters, method = c("IN", "BK", "PH")){
  n_fixed <- fixed(parameters)[1]

  if(method == "IN"){
    c2_cp <- function(z, n1){
      w <- rep(0,2)
      w[1] <- sqrt(n1 / n_fixed)
      w[2] <- sqrt( 1 - w[1]^2 )
      p <- (qnorm(1 - parameters$alpha) - w[1]*z)/w[2]
      return(p)
    }
  }

  if(method == "BK"){
    c_alpha <- function(parameters){
      exp(-0.5 * qchisq(1 - parameters$alpha, df=4))
    }

    c2_cp <- function(z, n1){
      qnorm(1 - c_alpha(parameters) / pnorm(-z) )
    }

  }

  n2_cp <- function(z, n1, cond_power){
    # Dringend verbessern!! S. Wassmer Buch Tabelle 2.5 und S. 104
    A <- function(z){
      u <- ifelse(parameters$alpha == 0.025, 2.168, 1.855)
      return(1 - pnorm(sqrt(2) * u - z))
    }
    p <- 2* ( (qnorm(A(z)) + qnorm(1 - cond_power)) / parameters$mu )^2
    # p <- n1 * ( (qnorm(A(z)) + qnorm(1 - cond_power)) / z )^2
    return(p)
  }


  score_cp <- function(n1, cf, ce){
    h = (ce - cf) / 10
    ww = seq(cf, ce, h)
    alph = c(1, rep(2, 9), 1)

    # z <- (z, alpha)
    f <- function(z, al){
      q <- al * n2_cp(z, n1, 1-parameters$beta) * dnorm(z - sqrt(n1) * parameters$mu)
      return(q)
    }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <- (h/2) * sum(y)
    p <- p + n1
    return(p)
  }

  power_lack <- function(n1, cf, ce){
    h = (ce - cf) / 10
    ww = seq(cf, ce, h)
    alph = c(1, rep(2, 9), 1)

    # Define type two error
    f <- function(z, al){
      al * pnorm(c2_cp(z, n1) - sqrt(n2_cp(z, n1, 1-parameters$beta)) * parameters$mu) * dnorm(z - sqrt(n1) * parameters$mu)
    }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <- parameters$beta - (h/2) * sum(y) - pnorm(cf - sqrt(n1) * parameters$mu)
    return(p)
  }

  type_one <- function(n1, cf, ce){
    h = (ce - cf) / 10
    ww = seq(cf, ce, h)
    alph = c(1, rep(2, 9), 1)

    # Define type one error
    f <- function(z, al){
      al * pnorm(c2_cp(z, n1)) * dnorm(z)
    }
    y <- apply(cbind(ww, alph), 1, function(x) f(x[1], x[2]))
    p <- 1 - (h/2) * sum(y) - pnorm(cf)
    return(p)
  }


  optimum <- nloptr::nloptr(
     x0          = c(n_fixed/2, 0, 2),
     eval_f      = function(x) score_cp(x[1], x[2], x[3]),
     eval_g_ineq = function(x) c( -power_lack(x[1], x[2], x[3]),
                                  type_one(x[1], x[2], x[3]) - parameters$alpha),
     lb = c(0.1 * n_fixed, 0, qnorm(1-parameters$alpha)),
     ub = c(0.9 * n_fixed, qnorm(1-parameters$alpha)-0.2, 4),
     opts = list(
       algorithm = "NLOPT_LN_COBYLA",
       xtol_rel = 0.0001,
       maxeval = 999999,
       maxtime = 16200
     )
   )

   n1_out <- optimum$solution[1]
   cf_out <- optimum$solution[2]
   ce_out <- optimum$solution[3]


  # noch machen!

  n2_out <- function(z){
    p <-ifelse(z <= ce_out && z >= cf_out, n2_cp(z, n1_out, 1-parameters$beta), 0)
    return(p)
  }
  c2_out <- function(z){
    if(z<= ce_out && z >= cf_out){
      p <- c2_cp(z, n1_out)
    } else if( z >= ce_out){
      p <- Inf
    }else{
      p <- -Inf
    }
  }


  d <- design(cf_out, ce_out, n1_out, n2_out, c2_out)
  return(d)


}
