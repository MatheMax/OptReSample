#' Find the optimal design when using a t-approximation
#'
#' \code{t_design} computes the optimal design when a t-approximation is used.
#'
#' @param parameters Parameters specifying the design
#'
#' @return An object of class \link{design}
#'
#' @export

t_design <- function(parameters){

  # Define grid
  min <- -0.5
  max <- 3
  dis <- 0.25

  w = seq(min, max, dis)


   # Start values
   start <- rep( 0 , 2*length(w)+3 )
   start[length(w)+3] <- (fixed(parameters)[1]) / 2
   start[2] <- sqrt(start[length(w)+3]) * parameters$mu
   start[3:(length(w)+2)] <- seq( 2.5, 0, -2.5/(length(w)-1) )
   start[(length(w)+4) : length(start)] <- seq( 60, 10, -50/(length(w)-1) )

    score_min <- function(cf, ce, c2, n1, n2){ t_score(parameters, cf, ce, c2, n1, n2, w) }
    t_1 <- function(cf, ce, c2, n1, n2){ t_type_one(parameters, cf, ce, c2, n1, n2, w) }
    t_2 <- function(cf, ce, c2, n1, n2){ t_type_two(parameters, cf, ce, c2, n1, n2, w) }

  optimum <- nloptr::nloptr(
    x0          = start,
    eval_f      = function(x) {score_min(x[1], x[2], x[3 : (length(w)+2)], x[length(w)+3], x[(length(w)+4) : length(start)])},
    eval_g_ineq = function(x) { return( c( x[1] + 0.01 - x[2] ,
                                  t_1(x[1], x[2], x[3 : (length(w)+2)], x[length(w)+3], x[(length(w)+4) : length(start)]) -
                                    parameters$alpha ,
                                  t_2(x[1], x[2], x[3 : (length(w)+2)], x[length(w)+3], x[(length(w)+4) : length(start)]) -
                                    parameters$beta )) },
    lb = c( rep(min,length(w)+2) , rep(2,length(w)+1) ),
    ub = c( rep(max,length(w)+2) , rep(400,length(w)+1) ),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 0.00001,
      maxeval = 999999,
      maxtime = 16200
    )
  )
  cf <- optimum$solution[1]
  ce <- optimum$solution[2]
  c2 <- optimum$solution[3:(length(w)+2)]
  n1 <- optimum$solution[length(w)+3]
  n2 <- optimum$solution[(length(w)+4) : length(start)]

  n2_out <- function(z){
    # Define bound
    a <- seq(min,cf,dis)
    b <- seq(ce,max,dis)
    l <- length(w)-length(a)-length(b)+1
    x <- rep(0,l)
    v <- rep(0,l)
    for(i in 1:l){
      x[i] <- w[i+length(a)]
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
    l <- length(w)-length(a)-length(b)+1
    x <- rep(0,l)
    v <- rep(0,l)
    for(i in 1:l){
      x[i] <- w[i+length(a)]
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


