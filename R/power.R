#' Computing the conditional power for differnt true effect sizes
#'
#' \code{real_power} gives the power of a design at a specified effect size
#'
#' In \code{parameters} one has to specify type I-, type II-error and effect size.
#'
#' @param effect Effect size on which one wants to calculate the power.
#' @param d An object of class \code{design}.
#' @param parameters Parameters speficying the situation. See \link{parameters} for details.

real_power<-function(effect,d,parameters){
  mualt=effect #the true effect size
  mu0=parameters$mu0
  sigma=parameters$sigma

  n1 <- d$n1
  cf <- d$cf
  ce <- d$ce

  # Define nods
  N=6
  h=(ce-cf)/(4*N)
  x=nodes(cf,ce,N)
  omega=rep(0,4*N+1)
  omega[1]=7
  w=c(32,12,32,14)
  omega[-1]=rep(w,N)
  omega[4*N+1]=7
  y=rep(0,4*N+1)
  for(i in 1:(4*N+1)){
    n <- as.numeric(d$n2(x[i]))
    c <- as.numeric(d$c2(x[i]))
    y[i] = pnorm( c - sqrt(abs(n)) * (mualt-mu0)/sigma ) * dnorm( x[i] - sqrt(abs(n1)) * (mualt-mu0)/sigma )
  }
  p <- (2*h)/45*(t(omega)%*%y)
  f <- 1 - pnorm ( cf - sqrt(abs(n1)) * (mualt-mu0)/sigma ) - p
  return(f)
}

#' Plot power dependent on true effect
#'
#' \code{plot_power} plots the power curve for the true effect.
#'
#' The true effect size has to be specified in \code{parameters} as well as type I - and type II - error constraints.
#'
#' @param d An object of class \code{design}.
#' @param parameters Parameters speficying the situation. See \link{parameters} for details.


plot_power<-function(d,parameters){
  mu0=parameters$mu0
  nu=parameters$mualt
  dis=(nu-mu0)/20
  z = seq(mu0-5*dis,1.5*nu,dis)
  y = rep(0,length(z))
  for(i in 1:length(z)){
   y[i] <- real_power(z[i],d,parameters)
  }
  out <- data.frame(data.matrix(cbind(z,y)))
  names(out)<-c("true effect","power")
  ggplot2::ggplot(out,ggplot2::aes(z, y)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1-parameters$beta, color = "red") +
    ggplot2::geom_vline(xintercept = parameters$mualt, color = "blue") +
    ggplot2::scale_x_continuous("true effect", breaks = seq(0, 10, .1)) +
    ggplot2::labs(title="Power",x="true effect",y="power")+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )

}


#' Compute the power of a design
#'
#' \code{opt_power} gives the power of a design on the desired point alternative.
#'
#' It equals \link{real_power} when \code{effect = parameters$mualt}.
#'
#' @param d An object of class \code{design}.
#' @param parameters Parameters speficying the situation. See \link{parameters} for details.

opt_power <- function(d,parameters) {
  real_power(parameters$mualt,d,parameters)
}
