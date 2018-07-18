#' Compute the design by Jennison and Turnbull (2015)
#'
#' \code{jt_design} computes the design proposed by Jennison and Turnbull (2015).
#' It is based on the inverse normal combination test.
#' The n_2-function is termined in order to minimize expected sample size under the alternative.
#'
#' Note that the sample sizes must be given per group. The total sample size is 2*(n_1 + n_2).
#'
#' @param theta Effect size for which the expected sample size should be minimized
#' @param sigma Known standard deviation
#' @param n1 First stage sample size. Write \code{NA} to obtain the optimal n_1.
#' @param lb_n2 Lower bound for the stage two sample size due to e.g. delayed response
#' @param alpha Maximal type I error rate
#' @param beta Maximal type II error rate
#'
#' @export


jt_design <- function(theta, sigma, n1=NA, lb_n2, alpha, beta){
  # Standardize effect
  th <- theta / sqrt(2) / sigma
  n_fixed <- fixed(list(alpha=alpha, beta=beta, mu=th))[1]
  x <- seq(-1, 4, 0.01)

  # lb_n2 > 0, n1 predefined
    if(is.numeric(n1)){
      # Combination test c2-function
      w_1 <- sqrt( n1 / 221 )
      w_2 <- sqrt( 1 - w_1^2 )
      c2_jt <- function(z){ (qnorm(1 - alpha) - w_1 * z) / w_2 }

      # n_2 function by JT
      n_2_test <- function(z, lambda){
        f <- function(n){
          log(8 * pi / (lambda * th)^2) + log(n) + (c2_jt(z) - sqrt(n) * th)^2
        }
        res <- try(uniroot(f, c(max(lb_n2, 1), 2000))$root, silent=T)
        if(class(res)=="try-error"){ res = lb_n2 }
        return(as.numeric(res))
      }

      # Find  lambda
      # Define power
      pow <- function(lambda){
        cp <- function(z){ 1 - pnorm( c2_jt(z) - sqrt(n_2_test(z, lambda)) * th ) }
        g <- function(z){ cp(z) * dnorm(z - sqrt(n1) * th) }
        p <- integrate(Vectorize(g), -Inf, Inf)$value
        return(p)
      }

    #lambda_opt <- 8 * sigma^2
    #while(pow(lambda_opt) < 1 - beta){lambda_opt <- lambda_opt + 1}
    lambda_opt <- 463

     # Define design
     n2 <- rep(0,length(x))
     for(i in 1:length(x)){
       n2[i] <- n_2_test(x[i], lambda_opt)
     }

     n1_out <- n1
     c2_out <- c2_jt

    } else{
    # lb_n2 > 0 , n1 optimized
      # Combination test c2-function
      c2_jt <- function(z, n1){
        w_1 <- sqrt( n1 / 221 )
        w_2 <- sqrt( 1 - w_1^2 )
        return((qnorm(1 - alpha) - w_1 * z) / w_2 )
      }

      # n_2 function by JT
      n_2_test <- function(z, n1, lambda){
        f <- function(n){
          log(8 * pi / (lambda * th)^2) + log(n) + (c2_jt(z, n1) - sqrt(n) * th)^2
        }
        res <- try(uniroot(f, c(lb_n2, 2000))$root, silent=T)
        if(class(res)=="try-error"){ res = lb_n2 }
        return(as.numeric(res))
      }

      # Find optimal lambda
      # Define objective criterion
      oc <- function(n1, lambda){
        f <- function(z){ (n_2_test(z, n1, lambda) - lb_n2) * dnorm(z - sqrt(n1) * th) }
        p <- integrate(Vectorize(f), -Inf, Inf)$value + n1 + lb_n2
        return(p)
      }
      # Define power
      pow <- function(n1, lambda){
        cp <- function(z){ 1 - pnorm( c2_jt(z, n1) - sqrt(n_2_test(z, n1, lambda)) * th ) }
        g <- function(z){ cp(z) * dnorm(z - sqrt(n1) * th) }
        p <- integrate(Vectorize(g), -Inf, Inf)$value
        return(p)
      }
'
      optimum <- nloptr::nloptr(
        x0          = c(n_fixed/2, max(8 * sigma^2, 400)),
        #x0 = c(100, 400),
        eval_f      = function(x) oc(x[1], x[2]),
        eval_g_ineq = function(x) {1 - beta - pow(x[1], x[2])},
        lb = c(1, 10),
        ub = c(floor(n_fixed), Inf),
        opts = list(
          algorithm = "NLOPT_LN_COBYLA",
          xtol_rel = 0.0001,
          maxeval = 999999,
          maxtime = 16200
        )
      )
      n1_out <- optimum$solution[1]
      lambda_opt <- optimum$solution[2]
'
      lambda_opt <- 8 *sigma^2

      optimum <- nloptr::nloptr(
        x0          = 90,
        #x0 = c(100, 400),
        eval_f      = function(x) oc(x, lambda_opt),
        eval_g_ineq = function(x) {1 - beta - pow(x, lambda_opt)},
        lb = 1,
        ub = floor(n_fixed),
        opts = list(
          algorithm = "NLOPT_LN_COBYLA",
          xtol_rel = 0.0001,
          maxeval = 999999,
          maxtime = 16200
        )
      )
      n1_out <- optimum$solution[1]

      # Define c2
      c2_out <- function(z){ c2_jt(z, n1_out) }

      # Define n2
      n2 <- rep(0,length(x))
      for(i in 1:length(x)){
        n2[i] <- n_2_test(x[i], n1_out, lambda_opt)
      }
    }

    # Define design
    cf <-  max( 0.5 * (x[min(which(n2 > lb_n2))] + x[min(which(n2 > lb_n2)) - 1]) , -1)
    ce <-  min( 0.5 * (x[max(which(n2 > lb_n2))] + x[max(which(n2 > lb_n2)) + 1]) , 3)
    n2_jt <- n2[which(n2 > lb_n2)]

    dis <- (round(ce, 1) - round(cf, 1)) / (length(n2_jt) - 1)
    u <- seq(round(cf, 1), round(ce, 1), dis)

    n2_out <- function(z){
      spl <- splinefun(u,n2_jt)
      p <- ifelse(cf <= z && z <= ce, spl(z), floor(lb_n2))
      return(p)
    }

  #}

  d <- design(cf, ce, n1_out, n2_out, c2_out)
  return(d)
}
