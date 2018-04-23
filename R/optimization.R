#' Compute an optimal adaptive design
#'
#' \code{optimal_design} returns an optimal two-stage adaptive design.
#'
#' The optimality criterion is given as the expected sample size under the null hypothesis.
#' If \code{standardized=T}, the values of \code{effect_null} and {sd} are set to their default, even if they were specified.
#'
#' The use of a t-approximation (\code{t_approx=T}) or a Lagrangian approach (\code{lagrange=T}) yields more precise results,
#' but makes the calculation slowlier.
#'
#' Note that there is no Lagrangian implementation of the t_approxmation. Hence, the value of \code{lagrange} will be
#' ignored, if \code{t_approx} is set to \code{T}.
#'
#' @param effect The effect on that the power is computed
#' @param alpha The maximal type I error rate
#' @param pow The minimal power
#' @param effect_null The effect under the null hypothesis. Default is 0.
#' @param sd The standard deviation. Default is 1.
#' @param standardized Logical. Is \code{effect} already standardized? Default is \code{T}.
#' @param t_approx Logical. Should a t-approximation be used? Default is \code{F}.
#' @param lagrange Logical. Should a Lagrangian procedure be used? Default is \code{T}.
#'
#' @return An object of class \link{design}


optimal_design <- function(effect, alpha, pow, effect_null=0, sd=1, standardized=T, t_approx=F, lagrange=T){
  beta <- 1-pow
  mu0 <- effect_null
  sigma <- sd

  if( standardized==T ){
    mu0 <- 0
    sigma <- 1
  }


  if( t_approx==T ){
    mu0 <- 0
    sigma <- 1
  }


  pars<-list(
    mualt=effect,
    mu0=mu0,
    sigma=sigma,
    alpha=alpha,
    beta=beta
  )

  if(t_approx == T){
    d <- t_design(pars)
  } else{
    if(lagrange == F){
      d <- direct_design(pars)
    } else{
      l <- find_lambda(pars)
      d <- lagrange_design(pars,l[1],l[2])
    }
  }

  return(d)

}
