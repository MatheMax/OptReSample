#' Type one error when using t-approximation
#'
#' \code{t_type_one} gives the type I error when a t-approximation is used.
#'
#' Mainly needed for \link{t_design}
#'
#' @param parameters Parameters specifying the design
#' @param cf Boundary for stopping for futility
#' @param ce Boundary for stopping for efficacy
#' @param n1 Stage one sample size
#' @param n2 n_2 values on the nodes
#' @param w nodes
#'
#' @export

t_type_one <- function(parameters, cf, ce, c2, n1, n2, w){
  f <- splinefun(w, c2)
  g <- splinefun(w, n2)
  N=12
  h = (ce - cf) / N
  x = seq(cf,ce,h)
  alpha=c(1, rep(2, (N-1)), 1)

  # x_a = c(x,alpha)
  to <- function(x_a){
    x_a[2] * pt( f(x_a[1]), df = max(g(x_a[1])-1,1) , ncp=0 ) * dt( x_a[1], df=n1-1 , ncp=0 )
  }


  y <- apply(cbind(x,alpha), 1, to)
  p <- (h/2) * sum(y)
  p <- 1 - pt( cf, df=n1-1 , ncp=0 ) - p
  return(p)
}
