#' Computing optimal c_2-values
#'
#' \code{c2} returns the optimal value of the c_2-function in the Lagrangian framework
#'
#'   @param parameters Parameters specifying the design
#'   @param z Z_1-value of the first stage
#'   @param n1 First stage sample size
#'   @param lambda1 Penalization parameter for type I error
#'   @param lambda2 Penalization parameter for type II error
#'   @param cf Boundary for stopping for futility
#'   @param ce Boundary for stopping for efficacy

c2 <- function(parameters, z, n1, lambda1, lambda2, cf, ce){
  n2 <- response(parameters, n1, lambda1, lambda2, z)
  q <- ( parameters$mu^2 * n2 - b(parameters, z, n1, lambda1, lambda2) ) / ( 2 * parameters$mu * sqrt(n2) )

  if(z<cf) {
    q=Inf
  }
  if(z>ce) {
    q=-Inf
  }

  return(q)
}


#'Plot the c2-function
#'
#' \code{plot_c2} plots the c_2 function of a design in dependence of the first stage z-value.
#'
#' @param d An object of class \link{design}

plot_c2 <- function(d) {
  dis = d$ce - d$cf
  h = dis / 30
  z = seq(d$cf,d$ce,h)
  y = rep(0,length(z))
  for(i in 1:length(z)){
    y[i] <- d$c2(z[i])
  }
  out <- data.frame(data.matrix(cbind(z,y)))
  names(out)<-c("z_1","c_2")
  ggplot2::ggplot(out, ggplot2::aes(z, y)) +
    ggplot2::geom_line() +
    ggplot2::labs(title="c_2(z_1)",x="z_1",y="c_2")+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}
