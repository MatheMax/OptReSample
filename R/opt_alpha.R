#' Compute the type I error rate of a design
#'
#' \code{opt_alpha} gives the type I error rate of a design.
#'
#' @param d An object of class \code{design}.
#' @export

opt_alpha <- function(d){
  n1 <- d$n1
  cf <- d$cf
  ce <- d$ce
  N=10
  h=(ce-cf)/(4*N)
  x=nodes(cf,ce,N)
  omega=rep(0,4*N+1)
  omega[1]=7
  w=c(32,12,32,14)
  omega[-1]=rep(w,N)
  omega[4*N+1]=7
  y=rep(0,4*N+1)
  for(i in 1:(4*N+1)){
    c <- d$c2(x[i])
    y[i] = pnorm(c) * dnorm(x[i])
  }
  p <- (2*h)/45*(t(omega)%*%y)
  f <- 1 - pnorm(cf) - p
  return(f)
}
