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

t_score <- function(parameters, cf, ce, c2, n1, n2, w){
  g <- splinefun(w, n2)

  N=12

  h = (ce - cf) / N
  x = seq(cf,ce,h)
  alpha=c(1,rep(2,(N-1)),1)

  y=rep(0,(N+1))
  for(i in 1:(N+1)){
    alpha[i] = alpha[i] * g(x[i])
    y[i] <- dt( x[i], df=n1-1, ncp = sqrt(n1) * parameters$mu )
  }
  p <- (h/2)*(t(alpha)%*%y)
  p <- p + n1
  return(p)
}

