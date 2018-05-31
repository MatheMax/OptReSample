#' Compute concrete n_2-value
#'
#' \code{response} computes the concrete n_2-value for a given first stage z_1-value.
#'
#'  Here, the Lagrangian framework is used du to the fact that the computation of n_2 can be reduced to a node-finding problem.
#'
#'  @param parameters Parameters specifying the design
#'  @param n1 First stage sample size
#'  @param lambda1 Penalization parameter for type I error
#'  @param lambda2 Penalization parameter for type II error
#'  @param z Z_1-value of the first stage
#'  @export


response <- function( parameters, n1, lambda1, lambda2, z ) {
  f <- function(v) {
    q <- 8*log( sqrt(8*pi) / ( abs(lambda2) * parameters$mu ) ) + 2 * b( parameters, z, n1, lambda1, lambda2 )
    q <- q + 4 * log( abs(v) ) + parameters$mu^2 * v + ( b( parameters, z, n1, lambda1, lambda2 ) / parameters$mu )^2 * (1/v)
    return(q)
  }
  min <- as.numeric( (-2 + sqrt( 4 + b(parameters,z,n1,lambda1,lambda2)^2 ) ) / parameters$mu^2 )
  w <- ifelse(is.finite(f(min)),f(min),1)
  if( w > 0 ){
    p <- 9999*lambda1*lambda2
  } else{
    # The following is actually sensible, but destroys the results...
    #if(f(1)>0 && min>1){
    #  p <- uniroot(f,c(1,min))$root
    #} else{
    ub <- min + 1
    fub <- f(ub)
    while(fub<0){
      ub <- ub + 1
      fub <- f(ub)
    }
    p <- uniroot(f,c(min,ub))$root
  }
  return(p)
}


#' Compute n_2 function
#'
#' Computes the optimal n_2-value for a given first stage z_1-value in the Lagrangian framework.
#'
#' \code{n2} differs from \link{response} only so far that the values are set to 0 outside the interval (cf,ce).
#'
#'  @param parameters Parameters specifying the design
#'  @param z Z_1-value of the first stage
#'  @param n1 First stage sample size
#'  @param lambda1 Penalization parameter for type I error
#'  @param lambda2 Penalization parameter for type II error
#'  @param cf Boundary for stopping for futility
#'  @param ce Boundary for stopping for efficacy
#' @export


n2 <- function(parameters, z, n1, lambda1, lambda2, cf, ce ) {
  p <- 0
  if( z >= cf && z <= ce){
    p <- response(parameters, n1, lambda1, lambda2, z)
  }
  return(p)
}


#' Plot the n2 function
#'
#' \code{plot_n2} plots the n2-function of a design.
#'
#' @param d An object of class \link{design}.
#' @export

plot_n2 <- function(d) {
  h = (d$ce - d$cf) / 30
  z = seq(d$cf,d$ce,h)
  y = rep(0,length(z))
  for(i in 1:length(z)){
    y[i] <- d$n2(z[i])
  }
  out <- data.frame(data.matrix(cbind(z, y)))
  names(out)<-c('z_1', "n_2")
  ggplot2::ggplot(out,ggplot2::aes(z, y)) +
    ggplot2::geom_line() +
    ggplot2::labs(title=expression(n[2](z[1])), x=expression(z[1]), y=expression(n[2]))+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}
