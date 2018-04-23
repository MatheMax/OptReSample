#' Compute the optimal group sequential design
#'
#' \code{optimal_gsd} computes the group sequential design that is optimal w.r.t. expected sample size under the
#' alternative hypothesis and holdes presepcified type I and type II error constraints.
#'
#' @param parameters Parameters specifying the design


optimal_gsd <- function(parameters) {
  # expected sample size
  exp_n_gs<-function(parameters,cf,ce,n1,c2,n2){
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
      y[i] = n2 * dnorm( ( x[i] - sqrt(abs(n1))*a(parameters) ) )
    }
    p <- (2*h)/45*(t(omega)%*%y)
    f <- n1 + p
    return(f)
  }

  # probability to reject
  reject <- function(parameters,cf,ce,n1,c2,n2,mu){
    p <- pnorm ( ce - sqrt(abs(n1)) * (mu - parameters$mu0) / parameters$sigma )
    p <- p - pnorm ( cf - sqrt(abs(n1)) * (mu - parameters$mu0) / parameters$sigma )
    p <- p * pnorm ( c2 - sqrt(abs(n2)) * (mu - parameters$mu0) / parameters$sigma )
    p <- -p + 1 - pnorm( cf - sqrt(abs(n1)) * (mu - parameters$mu0) / parameters$sigma )
    return(p)
  }


  #Type I error
  type_one <- function(cf,ce,n1,c2,n2){
    reject(parameters,cf,ce,n1,c2,n2,parameters$mu0)
  }

  #Power
  power_gsd <- function(cf,ce,n1,c2,n2){
    reject(parameters,cf,ce,n1,c2,n2,parameters$mualt)
  }

  exp_n_opt <- function(cf,ce,n1,c2,n2){
    exp_n_gs(parameters,cf,ce,n1,c2,n2)
  }


  optimum <- nloptr::nloptr(
    x0          = c(1,2,30,1,30),
    eval_f      = function(x) exp_n_opt(x[1],x[2],x[3],x[4],x[5]),
    eval_g_ineq = function(x) c( x[1] + 0.01 - x[2] ,type_one(x[1],x[2],x[3],x[4],x[5])-parameters$alpha,
                                 1 - parameters$beta - power_gsd(x[1],x[2],x[3],x[4],x[5])),
    lb = c(-3,-3,1,-3,1),
    ub = c(5,5,1000,5,1000),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.00001,
      maxeval = 3000
    )
  )
  cf <- optimum$solution[1]
  ce <- optimum$solution[2]
  n1 <- optimum$solution[3]
  c2 <- optimum$solution[4]
  n2 <- optimum$solution[5]
  n2_out <- function(z){ifelse(cf <= z && z <= ce , n2 ,0)}
  c2_out <- function(z){
    p=0
    if(z<cf){p=Inf}
    if(z>ce){p=-Inf}
    if(cf <= z && z <= ce){p=c2}
    return(p)
  }

  d <- design(cf,ce,n1,n2_out,c2_out)
  return(d)
}
