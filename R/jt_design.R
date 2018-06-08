#' Compute the design by Jennison and Turnbull (2015)
#'
#' \code{jt_design} computes the design proposed by Jennison and Turnbull (2015).
#' It is based on the inverse normal combination test.
#' The n_2-function is termined in order to minimize expected sample size under the alternative.
#'
#' @param theta Effect size for which the expected sample size should be minimized
#' @param sigma Known standard deviation
#' @param n1 First stage sample size
#' @param lb_n2 Lower bound for the stage two sample size due to e.g. delayed response
#' @param alpha Maximal type I error rate
#' @param beta Maximal type II error rate
#'
#' Note that the sample sizes must be given per group. Hence, the total sample size is 2*(n_1 + n_2).
#'
#' @export


jt_design <- function(theta, sigma, n1, lb_n2, alpha, beta){
  # Standardize effect
  th <- theta / sqrt(2) / sigma

  # Combination test c2-function
  pars <- list(alpha=alpha, beta=beta, mu=th)
  n_fixed <- fixed(pars)[1]
  w_1 <- sqrt( n1 / n_fixed )
  w_2 <- sqrt( 1 - w_1^2 )
  c2_jt <- function(z){ (qnorm(1 - alpha) - w_1 * z) / w_2 }

  # n_2 function by JT
  n_2_test <- function(z, lambda){
    f <- function(n){
      log(8 * pi / (lambda * th)^2) + log(n) + (c2_jt(z) - sqrt(n) * th)^2
    }
    res <- try(uniroot(f, c(lb_n2, 2000))$root, silent=T)
    if(class(res)=="try-error"){ res = lb_n2 }
    return(as.numeric(res))
  }

  # Find optimal lambda
    # Define objective criterion
    oc <- function(lambda){
      f <- function(z){ (n_2_test(z,lambda) - lb_n2) * dnorm(z - sqrt(n1) * th) }
      p <- integrate(Vectorize(f), -Inf, Inf)$value
      return(p)
    }
    # Define power
    pow <- function(lambda){
      cp <- function(z){ 1 - pnorm( c2_jt(z) - sqrt(n_2_test(z, lambda)) * th ) }
      g <- function(z){ cp(z) * dnorm(z - sqrt(n1) * th) }
      p <- integrate(Vectorize(g), -Inf, Inf)$value
      return(p)
    }

  optimum <- nloptr::nloptr(
      x0          = 8 * sigma^2,
      eval_f      = function(x) oc(x),
      eval_g_ineq = function(x) {1 - beta - pow(x)},
      lb = 0.1,
      ub = Inf,
      opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 0.0001,
        maxeval = 999999,
        maxtime = 16200
      )
  )
  lambda_opt <- optimum$solution

  # Define design
  x <- seq(-1, 3, 0.1)


  n2 <- rep(0,length(x))
  for(i in 1:length(x)){
    n2[i] <- n_2_test(x[i], lambda_opt)
  }


  cf <-  max( 0.5 * (x[min(which(n2 > lb_n2))] + x[min(which(n2 > lb_n2)) - 1]) , -1)
  ce <-  min( 0.5 * (x[max(which(n2 > lb_n2))] + x[max(which(n2 > lb_n2)) + 1]) , 3)
  n2_out <- n2[which(n2 > lb_n2)]

  dis <- (round(ce, 1) - round(cf, 1)) / (length(n2_out) - 1)
  u <- seq(round(cf, 1), round(ce, 1), dis)

  n2_jt <- function(z){
    spl <- splinefun(u,n2_out)
    p <- ifelse(cf <= z && z <= ce, spl(z), lb_n2)
    return(p)
  }


  d <- design(cf, ce, n1, n2_jt, c2_jt)
  return(d)
}
