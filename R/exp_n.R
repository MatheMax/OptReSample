#' Computing the expected sample size for an effect size
#'
#' \code{exp_n} computes the expected sample size for a design under a certain effect size.
#'
#' @param effect Effect size for which the expected sample size should be calculated
#' @param d An object of class \link{design}
#' @param parameters Parameters specifying the design
#' @export

exp_n<-function(effect, d, parameters){
  mualt=effect #the true standardized effect size

  n1 <- d$n1
  cf <- d$cf
  ce <- d$ce
  N=6
  h=(ce-cf)/(4*N)
  x=nodes(cf, ce, N)
  omega=rep(0,4*N+1)
  omega[1]=7
  w=c(32,12,32,14)
  omega[-1]=rep(w,N)
  omega[4*N+1]=7
  y=rep(0,4*N+1)
  for(i in 1:(4*N+1)){
    n <- d$n2(x[i])
    y[i] = n * dnorm( x[i] - sqrt(n1) * mualt )
  }
  p <- (2*h)/45*(t(omega)%*%y)
  f <- n1 + p
  return(f)
}

#' Plotting expected sample size
#'
#' \code{plot_exp_n} plots the expected sample size for different effect sizes
#'
#' @param d An object of class \link{design}
#' @param parameters Parameters specifying the design
#' @export

plot_exp_n<-function(d,parameters){
  nu = parameters$mu
  z = seq(0, 1.5*nu, nu/20)
  y = rep(0,length(z))
  for(i in 1:length(z)){
    y[i] <- exp_n(z[i],d,parameters)
  }
  out <- data.frame(data.matrix(cbind(z,y)))
  names(out)<-c("true effect","expected sample size")
  ggplot2::ggplot(out, ggplot2::aes(z, y)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = parameters$mu, linetype = "dotted") +
    ggplot2::scale_x_continuous("true standardized effect", breaks = seq(0, 10, .1)) +
    ggplot2::labs(title="Expected sample size",x="true effect",y="sample size") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}


#' Plotting total sample size
#'
#' \code{plot_n} plots the total sample size of a design dependent on the z-value of the first stage.
#'
#' @param d An object of class \link{design}
#' @param parameters Parameters specifying the design
#' @export

plot_n <- function(d, parameters) {
  dis = d$ce - d$cf
  h = dis / 30
  z = seq(d$cf,d$ce,h)
  y = rep(0,length(z))
  for(i in 1:length(z)){
    y[i] <- d$n2( z[i] ) + d$n1
  }
  out <- data.frame(data.matrix(cbind(z,y)))
  names(out)<-c("effect estimator","n")
  ggplot2::ggplot(out, ggplot2::aes(z, y)) +
    ggplot2::geom_point() +
    #scale_x_continuous("effect estimator", breaks = seq(d$cf,d$ce, .2)) +
    ggplot2::labs(title="Total sample size",x="z",y="sample size")+
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
}

