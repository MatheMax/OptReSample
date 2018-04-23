#' Computing conditional power
#'
#' \code{cond_power} computes the conditional power for a given z_1-value.
#'
#' Note that this function assumes a Gaussian distribution and a t-approximation has not been implemented yet.
#'
#' @param z Z_1-value of the first stage
#' @param d An object of class \code{design}
#' @param parameters Parameters secifying the design
#'

cond_power<-function(z,d,parameters){
  f <- 0
  if(d$cf <= z && z <= d$ce){
    f <- 1- pnorm( d$c2(z) - sqrt(d$n2(z)) * ( parameters$mualt - parameters$mu0 ) / parameters$sigma )
  }
  return(f)
}


#' Plotting conditional power
#'
#' \code{cond_power} plots the conditional power curve of a design.
#'
#' Note that this function assumes a Gaussian distribution and a t-approximation has not been implemented yet.
#'
#' @param d An object of class \code{design}
#' @param parameters Parameters secifying the design


plot_cond_power<-function(d,parameters){
  N = 5
  z = cnodes(d$cf,d$ce,N)
  y = rep(0,length(z))
  for(i in 1:length(z)){
   y[i] <- cond_power(z[i],d,parameters)
  }
  out <- data.frame(data.matrix(cbind(z,y)))
  names(out)<-c("z_1","conditional power")
  ggplot(out,aes(z, y)) +
    geom_line() +
    geom_hline(yintercept = 1-parameters$beta, color = "red") +
   # scale_x_continuous("z_1", breaks = seq(d$cf, d$ce, .2)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )

}


