#' Version of the score when using t-approximation
#'
#' \code{t_score} gives the version of the score when the t_approximation is used.
#'
#' Mainly needed for \link{t_design}
#'
#' @param parameters Parameters specifying the design
#' @param cf Boundary for stopping for futility
#' @param ce Boundary for stopping for efficacy
#' @param c2 c_2 values on the nodes
#' @param n1 Stage one sample size
#' @param n2 n_2 values on the nodes
#' @param w nodes
#'
#' @export

t_score <- function(parameters, cf, ce, c2, n1, n2, w){
  g <- splinefun(w, n2)

  N=12

  h = (ce - cf) / N
  x = seq(cf, ce, h)
  alpha=c(1,rep(2,(N-1)),1)

  # x_a = c(x, alpha)
  sc <- function(x_a){
    x_a[2] * g(x_a[1]) * dt( x_a[1], df=n1-1, ncp = sqrt(n1) * parameters$mu )
  }

  y <- apply(cbind(x,alpha), 1, sc)

  p <- (h/2) * sum(y)
  p <- p + n1
  return(p)
}

