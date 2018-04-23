#' Compute smooth designs via the direct method
#'
#' \code{direct_design_smooth} is a different way to compute optimal designs via a spline approach.
#'
#' This function uses \link{stage_two} to find optimal stage-two functions and optimizes then over the stage one
#' parameters n1, cf and ce. It is a more suitable way than the one that is used in \link{direct_design} because one only
#' evaluates the stage two-functions inside the interval (cf,ce) and not on a full grid.
#' This function yields smoother results than \link{direct_design}, but also takes more time.
#'
#' @param parameters Parameters specifying the design
#'
#' @return An object of class design.

direct_design_smooth <- function (parameters){
  s_min <- function(cf, ce, n1){ score_smooth(parameters, cf, ce, n1) }

  k <- optimal_gsd(parameters)
  start_cf <- k$cf
  start_ce <- k$ce
  start_n1 <- k$n1

  low <- c(qnorm(parameters$beta), qnorm(1-parameters$alpha), 1)
  up <- c(qnorm(1-parameters$alpha)-0.1, 4, Inf)


  optimum <- nloptr::nloptr(
    x0          = c( start_cf, start_ce, start_n1 ),
    eval_f      = function(x) s_min(x[1], x[2], x[3]),
    lb = low,
    ub = up,
    opts = list(
      algorithm = "NLOPT_LN_NELDERMEAD",
      # NELDERMEAD failes sometimes. therefore, alternatively:
      # algorithm = "NLOPT_LN_SBPLX",
      maxtime = 9999,
      xtol_rel = 0.0001,
      maxeval = 99999
    )
  )

  cf <- optimum$solution[1]
  ce <- optimum$solution[2]
  n1 <- optimum$solution[3]

  s2 <- stage_two(parameters,cf,ce,n1)
  c2 <- s2[1:(length(s2)/2)]
  n2 <- s2[(length(s2)/2+1):length(s2)]
  dis <- (round(ce,1)-round(cf,1))/(length(n2)-1)
  u <- seq(round(cf,1), round(ce,1), dis)

  n2_out <- function(z){
    spl <- splinefun(u,n2)
    p <- ifelse(cf <= z && z <= ce, spl(z) ,0)
    return(p)
  }
  c2_out <- function(z){
    spl <- splinefun(u,c2)
    if(z<cf){ p=Inf }
    if(z>ce){ p=-Inf }
    if(cf <= z && z <= ce){
      p <- ifelse(cf <= z && z <= ce, spl(z) ,0)
    }
    return(p)
  }

  d <- design(cf,ce,n1,n2_out,c2_out)
  return(d)
}



