#' Creating design objects
#'
#' \code{design} creats an object of class \code{design}.
#'
#' @param cf Boundary for stopping for futility
#' @param ce Boundary for stopping for efficacy
#' @param n1 First stage sample size
#' @param n2 Function of z which gives the stage two sample size
#' @param c2 Function of z which gives the stage two rejection boundary
#'
#' @return An object of class \code{design}

design <- function(
  cf,
  ce,
  n1,
  n2,
  c2
) {
  if (!is.numeric(n1)) {
    stop("n1 must be integer")
  }
  return(
    structure(list(
      n1 = n1,
      cf = cf,
      ce = ce,
      n2  = n2,
      c2  = c2
    ), class = c("Design"))
  )
}


#' Creating parameters object
#'
#' \code{parameters} creats an object of class \code{parameters} which is needed for many functions in this package.
#'
#' If \code{mu0} and \code{sigma} are not speficied, \code{mualt} is assumed to be the standardized effect size
#' of the alternative.
#'
#' @param mualt Effect value of the alternative hypothesis
#' @param alpha Maximal type I error rate
#' @param beta Maximal type II error rate
#' @param mu0 Effect value of the null hypothesis. Default value is 0.
#' @param sigma Standard deviation. Default value is 1.
#'
#' @return An object of class \code{parameters.}

parameters <- function(
  mualt,
  alpha,
  beta,
  mu0 = 0,
  sigma = 1
) {
    return(
    structure(list(
      mualt = mualt,
      mu0 = mu0,
      sigma = sigma,
      alpha = alpha,
      beta = beta
    ), class = c("Parameters"))
  )
}
