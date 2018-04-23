#' Computing the expected sample size for an effect size
#'
#' \code{exp_n} computes the expected sample size for a design under a certain effect size.
#'
#' @param effect Effect size for which the expected sample size should be calculated
#' @param d An object of class \link{design}
#' @param parameters Parameters specifying the design

exp_n<-function(effect,d,parameters){
  mualt=effect #the true effect size
  mu0=parameters$mu0
  sigma=parameters$sigma

  n1 <- d$n1
  cf <- d$cf
  ce <- d$ce
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
    n <- d$n2(x[i])
    y[i] = n * dnorm( ( x[i] - sqrt(abs(n1))*(mualt-mu0)/sigma ) )
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

plot_exp_n<-function(d,parameters){
  mu0=parameters$mu0
  mualt=parameters$mualt
  dis=(mualt-mu0)
  z = seq(mu0-dis/2,mualt+dis,dis/20)
  y = rep(0,length(z))
  for(i in 1:length(z)){
    y[i] <- exp_n(z[i],d,parameters)
  }
  out <- data.frame(data.matrix(cbind(z,y)))
  names(out)<-c("true effect","expected sample size")
  ggplot(out,aes(z, y)) +
    geom_line() +
    geom_vline(xintercept = parameters$mualt, color = "red") +
    scale_x_continuous("true effect", breaks = seq(0, 10, .1)) +
    labs(title="Expected sample size",x="true effect",y="sample size") +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
}


#' Plotting total sample size
#'
#' \code{plot_n} plots the total sample size of a design dependent on the z-value of the first stage.
#'
#' @param d An object of class \link{design}
#' @param parameters Parameters specifying the design

plot_n <- function(d,parameters) {
  dis = d$ce - d$cf
  h = dis / 30
  z = seq(d$cf,d$ce,h)
  y = rep(0,length(z))
  for(i in 1:length(z)){
    y[i] <- d$n2( z[i] ) + d$n1
  }
  out <- data.frame(data.matrix(cbind(z,y)))
  names(out)<-c("effect estimator","n")
  ggplot(out,aes(z, y)) +
    geom_point() +
    #scale_x_continuous("effect estimator", breaks = seq(d$cf,d$ce, .2)) +
    labs(title="Total sample size",x="z",y="sample size")+
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
}

