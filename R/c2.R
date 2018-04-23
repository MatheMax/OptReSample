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

c2 <- function(parameters,z,n1,lambda1,lambda2,cf,ce){
  N=20
  no <- nodes(cf,ce,N)
  resp <- rep(0,length(no))
  for(i in 1:length(no)){
    resp[i] <- response(parameters,n1,lambda1,lambda2,no[i])
  }
  p <- rep(0,length(no))
  for(i in 1:length(no)){
  p[i] <- ( a(parameters)^2 * resp[i] - b(parameters,z,n1,lambda1,lambda2) ) / ( 2 * a(parameters) * sqrt(abs(resp[i])) )
  }
  c <- splinefun(no,p)
  q = c(z)
  if(z<cf) {
    q=Inf
  }
  if(z>ce) {
    q=-Inf
  }
  if (is.na(q)) {
    stop("response must not be missing")
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
  ggplot2::ggplot(out,aes(z, y)) +
    geom_line() +
    labs(title="c_2(z_1)",x="z_1",y="c_2")+
    theme_bw() +
    theme(
      panel.grid = element_blank()
    )
}
