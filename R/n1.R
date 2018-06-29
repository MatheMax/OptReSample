#' Find the optimal stage one sample size
#'
#' \code{n1} computes the stage one sample size of an optimal adaptive design
#'
#' @param parameters The parameters (alpha, power, standardized effect) whith which you want to build your design
#' @param lambda1,lambda2 The Lagrange penalization parameters
#' @export


n1 <- function(parameters, lambda1, lambda2){
  s <- function(n){
    c <- c_early(parameters,n,lambda1,lambda2)
    p <- score(parameters, c[1], c[2], n, lambda1, lambda2)
    return(p)
  }

  opt <- nloptr::nloptr(
    x0        =  (fixed(parameters)[1]) / 2,
    eval_f  = function(x) { s(x) },
    lb = 1,
    ub = 1000,
    opts = list(
      algorithm = "NLOPT_LN_BOBYQA",
      xtol_rel = 0.0001,
      maxeval = 50000
    )
  )

  return(round(opt$solution[1]))
}
