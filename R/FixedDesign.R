#' Compute the fixed design
#'
#' \code{fixed} computes the fixed design.
#'
#' @param parameters Parameters specifying the design



fixed <- function(parameters){
  score_fixed <- function(parameters, n, c, lambda1, lambda2){
    p <- n + lambda1 * ( 1 - pnorm ( c ) - parameters$alpha )
    p <- p + lambda2 * ( - parameters$beta + pnorm( c - sqrt(n) * parameters$mu ) )
  }


  optimal <- nloptr::nloptr(
    x0          = c(40,2,120,40),
    eval_f      = function(x) score_fixed(parameters, x[1], x[2], x[3], x[4]),
    lb = c(1,-5,0,0),
    eval_g_ineq = function(x) c((1 - pnorm ( x[2] ) - parameters$alpha),
                                - parameters$beta + pnorm( x[2] - sqrt(x[1]) * parameters$mu )),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.00001,
      maxeval = 500
    )
  )

  n <- optimal$solution[1]
  c <- optimal$solution[2]
  lambda1 <- optimal$solution[3]
  lambda2 <- optimal$solution[4]

  return(c(n,c,lambda1,lambda2))
}
