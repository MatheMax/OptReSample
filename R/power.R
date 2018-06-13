#' Computing the conditional power for differnt true effect sizes
#'
#' \code{real_power} gives the power of a design at a specified effect size
#'
#' In \code{parameters} one has to specify type I-, type II-error and effect size.
#'
#' @param effect Effect size on which one wants to calculate the power.
#' @param d An object of class \code{design}.
#' @param parameters Parameters speficying the situation. See \link{parameters} for details.
#'
#' @export

real_power<-function(effect, d, parameters){
  # Define nodes
  N = 6
  h = (d$ce - d$cf) / (4 * N)
  x = nodes(d$cf, d$ce, N)
  omega = rep(0, 4 * N+1)
  omega[1] = 7
  w = c(32, 12, 32, 14)
  omega[-1] = rep(w, N)
  omega[4 * N+1] = 7

  # Compute Power
  y=rep(0, 4*N+1)
  for(i in 1:(4*N+1)){
    y[i] = pnorm( d$c2(x[i]) - sqrt(d$n2(x[i])) * effect ) * dnorm( x[i] - sqrt(d$n1) * effect )
  }
  p <- (2*h)/45*(t(omega)%*%y)
  f <- 1 - pnorm ( d$cf - sqrt(d$n1) * effect ) - p
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
#'
#' @export


plot_power<-function(d, parameters){
  nu = parameters$mu
  z = seq(0, 1.5*nu, nu/20)
  y = rep(0,length(z))
  for(i in 1:length(z)){
   y[i] <- real_power(z[i], d, parameters)
  }
  out <- data.frame(data.matrix(cbind(z, y)))
  names(out)<-c("true effect","power")
  ggplot2::ggplot(out,ggplot2::aes(z, y) ) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 1-parameters$beta, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = parameters$mu, linetype = "dotted") +
    ggplot2::scale_x_continuous("true standardized  effect", breaks = seq(0, 10, .1)) +
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
#'
#' @export

opt_power <- function(d, parameters) {
  real_power(parameters$mu, d, parameters)
}
