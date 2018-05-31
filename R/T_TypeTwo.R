#' Type two error when using t-approximation
#'
#' \code{t_type_two} gives the type II error when a t-approximation is used.
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

t_type_two <- function(parameters, cf, ce, c2, n1, n2, w){
  f <- splinefun(w,c2)
  g <- splinefun(w,n2)
  N=12
  h = (ce - cf) / N
  x = seq(cf,ce,h)
  alpha=c(1,rep(2,(N-1)),1)

  # x_a = c(x,alpha)
  tt <- function(x_a){
    x_a[2] *
      pt( f(x_a[1]), df = max(g(x_a[1])-1,1) , ncp = sqrt(g(x_a[1]))*parameters$mu ) *
      dt( x_a[1], df=n1-1, ncp = sqrt(n1)*parameters$mu )
  }

  y <- apply( cbind(x,alpha), 1, tt)
  p <- (h/2) * sum(y)
  p <- p + pt( cf, df=n1-1, ncp = sqrt(n1)*parameters$mu )
  return(p)
}
