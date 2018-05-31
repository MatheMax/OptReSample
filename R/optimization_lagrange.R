#' Find the optimal Lagrangian design
#'
#' Compute the optimal design for predefined Lagrangian penalty terms.
#'
#' Note that this design does not hold type I or type II error constraints necessarily.
#' To find values of lambda1 and lambda2 which yield the exact error constraints use \link{find_lambda}.
#'
#' @param parameters The parameters (alpha, power, etc.) you want to use for your design
#' @param lambda1 Penalization parameter for type I error
#' @param lambda2 Penalization parameter for type II error
#'
#' @return An object of class \code{design}
#'
#' @export


lagrange_design <- function(parameters, lambda1, lambda2) {
  n1 <- n1(parameters, lambda1, lambda2)
  c <- c_early(parameters, n1, lambda1, lambda2)
  cf <- c[1]
  ce <- c[2]
  n2_out <- function(z){ n2(parameters, z, n1, lambda1, lambda2, cf, ce) }
  c2_out <- function(z){ c2(parameters, z, n1, lambda1, lambda2, cf, ce) }

  d <- design(cf, ce, n1, n2_out, c2_out)
  return(d)
}



