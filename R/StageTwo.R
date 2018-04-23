#' Score version for smooth direct designs
#'
#' \code{score_direct_smooth} gives the version of the score that is needed for \link{stage_two}.
#'
#' @param parameters Parameters specifying the design
#' @param n1 First stage sample size
#' @param n2 n_2-values
#' @param h Distance between two nodes
#' @param N 4N+1 gives the number of nodes
#' @param w nodes inside the interval (cf,ce)

score_direct_smooth <- function(parameters,n1,n2,h,N,w){

  omega=rep(0,4*N+1)
  omega[1]=7
  se=c(32,12,32,14)
  omega[-1]=rep(se,N)
  omega[4*N+1]=7
  y=rep(0,4*N+1)
  for(i in 1:(4*N+1)){
    y[i] <- n2[i] * dnorm( w[i] - sqrt(n1)*a(parameters) )
  }
  p <- (2*h)/45*(t(omega)%*%y)


  p <- p + n1
  return(p)
}

#' Smooth score version
#'
#' \code{score_smooth} gives the version of the score that is needed for \link{score_direct_smooth}.
#'
#' @param Parameters Parameters specifying the design
#' @param cf Boundary for stopping for futility
#' @param ce Boundary for stopping for efficacy
#' @param n1 Stage one sample size

score_smooth <- function(parameters,cf,ce,n1){
  N=3
  h=(ce-cf)/(4*N)
  w <- seq(cf,ce,h)

  s2 <- stage_two(parameters,cf,ce,n1)
  n2 <- s2[ (length(s2)/2 + 1) : length(s2)]

  p <- score_direct_smooth(parameters,n1,n2,h,N,w)
  return(p)
}




#' Type I version for smooth direct designs
#'
#' \code{type_one_smooth} gives the version of the type I error that is needed for \link{stage_two}.
#'
#' @param parameters Parameters specifying the design
#' @param cf Boundary for stopping for futility
#' @param c2 c_2-values
#' @param h Distance between two nodes
#' @param N 4N+1 gives the number of nodes
#' @param w nodes inside the interval (cf,ce)


type_one_smooth <- function(parameters,cf,c2,h,N,w){

  omega=rep(0,4*N+1)
  omega[1]=7
  se=c(32,12,32,14)
  omega[-1]=rep(se,N)
  omega[4*N+1]=7
  y=rep(0,4*N+1)
  for(i in 1:(4*N+1)){
    y[i] <- pnorm(c2[i]) * dnorm(w[i])
  }
  p <- (2*h)/45*(t(omega)%*%y)


  p <- 1 - pnorm(cf) - p
  return(p)
}


#' Type II version for smooth direct designs
#'
#' \code{type_two_smooth} gives the version of the type II error that is needed for \link{stage_two}.
#'
#' @param parameters Parameters specifying the design
#' @param cf Boundary for stopping for futility
#' @param c2 c_2-values
#' @param n1 First stage sample size
#' @param n2 n_2-values
#' @param h Distance between two nodes
#' @param N 4N+1 gives the number of nodes
#' @param w nodes inside the interval (cf,ce)


type_two_smooth <- function(parameters,cf,c2,n1,n2,h,N,w){

  omega=rep(0,4*N+1)
  omega[1]=7
  se=c(32,12,32,14)
  omega[-1]=rep(se,N)
  omega[4*N+1]=7
  y=rep(0,4*N+1)
  for(i in 1:(4*N+1)){
    y[i] <- pnorm(c2[i] - sqrt(abs(n2[i])) * a(parameters) ) * dnorm(w[i] - sqrt(n1) * a(parameters))
  }
  p <- (2*h)/45*(t(omega)%*%y)

  p <- p + pnorm(cf - sqrt(n1) * a(parameters))
  return(p)
}



#' Compute optimal stage two values for given first stage
#'
#' \code{stage_two} computes the functions c_2 and n_2 that hold the error constraints and are optimal w.r.t.
#' expected sample size under the alternative for a prespecified first stage.
#'
#' \code{stage_two} is the base of \link{direct_design_smooth}.
#'
#' @param parameters Parameters specifying the design
#' @param cf Boundary for stopping for futility
#' @param ce Boundary for stopping for efficacy
#' @param n1 First stage sample size
#'
#' @return A vector. The first half give the c_2, and the second the n_2-values, on a equidistance grid inside the
#' interval (cf,ce).

stage_two <- function(parameters,cf,ce,n1){

  N=3
  h=(ce-cf)/(4*N)
  w <- seq(cf,ce,h)

  k <- optimal_gsd(parameters)
  start_n2 <- rep( ceiling( k$n2( k$cf + (k$ce-k$cf) / 2 ) ) , length(w) )
  start_c2 <- rep( k$c2( k$cf + (k$ce-k$cf) / 2 ) , length(w) )


  score_min <- function(n2){ score_direct_smooth(parameters,n1,n2,h,N,w) }
  t_1 <- function(c2){ type_one_smooth(parameters,cf,c2,h,N,w) }
  t_2 <- function(c2,n2){ type_two_smooth(parameters,cf,c2,n1,n2,h,N,w) }

  optimum <- nloptr::nloptr(
    x0          = c(start_c2,start_n2),
    eval_f      = function(x) score_min(x[(length(w)+1) : (2*length(w))]),
    eval_g_ineq = function(x) c( t_1(x[1:length(w)]) - parameters$alpha,
                                 t_2(x[1:length(w)],x[(length(w)+1) : (2*length(w))]) - parameters$beta),
    lb = c(rep(-1,length(w)),rep(1,length(w))),
    ub = c(rep(4,length(w)),rep(Inf,length(w))),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.0001,
      maxeval = 99920000
    )
  )

  c2<-optimum$solution[1:length(w)]
  n2 <- optimum$solution[(length(w)+1) : (2*length(w))]

  r1 <- t_1(c2) / parameters$alpha
  r2 <- t_2(c2,optimum$solution[(length(w)+1) : (2*length(w))]) / parameters$beta

  if(abs(1-r1)<0.05 && abs(1-r2)<0.05){ n2 <- optimum$solution[(length(w)+1) : (2*length(w))]}
  else{n2 <- rep(99999,length(w))}

  return(c(c2,n2))
}

