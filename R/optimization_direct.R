#' Find the optimal design via direct optimization
#'
#' \code{direct_design} comuptes the optimal design via a spline approach.
#'
#' The z-space is partionated in a grid and the optimal c_2- and n_2 values are computed on every node.
#' The function is truncated if z is not element of (c_f,c_e) and inside this interval, a spline function is fitted.
#'
#' @param parameters Parameters specifying the design.
#'
#' @return An object of class \link{design}.
#'


direct_design <- function(parameters){

  min <- -1.5
  max <- 4.5
  dis <- 0.5

  nodes = seq(min, max, dis)


   # Startwerte
   start_cf <- 0
   start_n1<- (fixed(parameters)[1])/2
   start_ce <- fixed(parameters)[2]
   start_c2<- seq(2.5,0,-2.5/(length(nodes)-1))
   start_n2 <- seq(60,10,-50/(length(nodes)-1))
   start <- c(start_cf, start_ce, start_c2, start_n1, start_n2)

   low <- c( min+0.1, qnorm(1-parameters$alpha) - dis, rep((min+0.1), length(nodes)), rep(1, length(nodes)+1) )
   up <- c( qnorm(1-parameters$alpha), max-0.1, rep((max-0.1), length(nodes)), rep(Inf, length(nodes)+1) )


    score_min <- function(cf, ce, c2, n1, n2){ score_direct(parameters, cf, ce, nodes, c2, n1, n2) }
    t_1 <- function(cf, ce, c2){ type_one(parameters, cf, ce, nodes, c2) }
    t_2 <- function(cf, ce, c2, n1, n2){ type_two(parameters, cf, ce, nodes, c2, n1, n2) }

  optimum <- nloptr::nloptr(
    x0          = start,
    eval_f      = function(x) score_min(x[1], x[2], x[3:(length(nodes)+2)], x[length(nodes)+3], x[(length(nodes)+4) : length(start)]),
    eval_g_ineq = function(x) { return( c( #x[1] + 0.01 - x[2] ,
                                  t_1(x[1], x[2], x[3:(length(nodes)+2)]) - parameters$alpha,
                                  t_2(x[1], x[2], x[3:(length(nodes)+2)], x[length(nodes)+3], x[(length(nodes)+4) : length(start)]) - parameters$beta )) },
    lb = low,
    ub = up,
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.0001,
      maxeval = 99999,
      maxtime = 16200
    )
  )
  cf <- optimum$solution[1]
  ce <- optimum$solution[2]
  c2 <- optimum$solution[3:(length(nodes)+2)]
  n1 <- optimum$solution[length(nodes)+3]
  n2 <- optimum$solution[(length(nodes)+4) : length(start)]

  n2_out <- function(z){
    # Define bound
    a <- seq(min,cf,dis)
    b <- seq(ce,max,dis)
    l <- length(nodes)-length(a)-length(b)+1
    x <- rep(0,l)
    v <- rep(0,l)
    for(i in 1:l){
      x[i] <- nodes[i+length(a)]
      v[i] <- n2[i+length(a)]
    }

    spl <- splinefun(x,v)
    p <- ifelse(cf <= z && z <= ce,spl(z),0)
    return(p)
  }

  c2_out <- function(z){
    # Define bound
    a <- seq(min,cf,dis)
    b <- seq(ce,max,dis)
    l <- length(nodes)-length(a)-length(b)+1
    x <- rep(0,l)
    v <- rep(0,l)
    for(i in 1:l){
      x[i] <- nodes[i+length(a)]
      v[i] <- c2[i+length(a)]
    }

    spl <- splinefun(x,v)
    if(z<cf){ p=Inf }
    if(z>ce){ p=-Inf }
    if(cf <= z && z <= ce){
      p <- spl(z)
    }
    return(p)
  }

  d <- design(cf,ce,n1,n2_out,c2_out)
  return(d)
}

