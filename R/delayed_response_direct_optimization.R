#' Compute designs for delayed response situation
#'
#' @export


dr_design <- function(theta, sigma, n1=NA, lb_n2, alpha, beta){
  # Define stuff for integration
  w <- seq(-1,3,0.2)
  N <- (length(w) - 1 ) /4
  h <- (max(w) - min(w)) / (4*N)
  omega = rep(0,4*N+1)
  omega[1] =7
  se =c(32,12,32,14)
  omega[-1] = rep(se,N)
  omega[4*N+1] = 7

  th = theta / sqrt(2) / sigma

  # Compute starting values
  gs <- optimal_gsd(list(alpha=alpha, beta=beta, mu=th))
  start_n2 <- rep(0, length(w))
  for(i in 1:length(w)){
    if(gs$cf <= w[i] && w[i] <= gs$ce){
      start_n2[i] <- ceiling( gs$n2( gs$cf + (gs$ce-gs$cf) / 2 ) ) - lb_n2
    }
  }
  start_c2 <- rep( gs$c2( gs$cf + (gs$ce - gs$cf) / 2 ) , length(w) )

  start<-c(start_c2, start_n2)

  # Score
  score_direct <- function(n1, n2){
    # x = (n2, w, omega)
    f <- function(x){
      x[3] * x[1] * dnorm(x[2] - sqrt(n1) * th)
    }
    y <- apply(cbind(n2, w, omega), 1 , f)
    p <- (2*h)/45 * sum(y)
    p <- p + n1 + lb_n2
    return(p)
  }

  # Type I error
  type_one <- function(c2){
    # x = (c2, w, omega)
    f <- function(x){
      x[3] * pnorm(x[1]) * dnorm(x[2])
    }
    y <- apply(cbind(c2, w, omega), 1 , f)

    p <- (2*h)/45 * sum(y)
    p <- 1 - pnorm(min(w)) - p
    return(p)
  }

  # Type II error
  type_two <- function(c2, n1, n2){
    # x = (c2, n2, w, omega)
    f <- function(x){
      x[4] * pnorm( x[1] - sqrt(lb_n2 + x[2]) * th ) * dnorm(x[3] - sqrt(n1) * th)
    }
    y <- apply(cbind(c2, n2, w, omega), 1 , f)

    p <- (2*h)/45 * sum(y)

    p <- p +  pnorm( min(w) - sqrt(n1) * th)
    return(p)
  }


  # Obtain the same power and type one error as JT
  jenturn <- jt_design(theta, sigma, 104, lb_n2, alpha, beta)

  c2_test <- rep(0,length(w))
  n2_test <- rep(0,length(w))

  for(i in 1:length(w)){
    c2_test[i] <- jenturn$c2(w[i])
    n2_test[i] <- jenturn$n2(w[i]) - lb_n2
  }

  al <- type_one(c2_test)
  be <- type_two(c2_test, 104, n2_test)



  # Compute second stage
  stage_two <- function(n1){

    score_min <- function(n2){ score_direct(n1,n2) }
    t_2 <- function(c2,n2){ type_two(c2,n1,n2) }

    optimum <- nloptr::nloptr(
      x0          = start,
      eval_f      = function(x) score_min(x[(length(w)+1) : (length(start))]),
      eval_g_ineq = function(x) c( type_one(x[1:length(w)]) - al,
                                   t_2(x[1:length(w)], x[(length(w)+1) : (length(start))]) - be),
      lb = c(rep(-1,length(w)), rep(0,length(w))),
      ub = c(rep(4,length(w)), rep(300,length(w))),
      opts = list(
        algorithm = "NLOPT_LN_COBYLA",
        xtol_rel = 0.0001,
        maxeval = 9999999
      )
    )

    c2 <- optimum$solution[1:length(w)]

    r1 <- type_one(c2) / al
    r2 <- t_2(c2,optimum$solution[(length(w)+1) : (length(start))]) / be

    if(abs(1-r1)<0.05 && abs(1-r2)<0.05){
      n2 <- round(optimum$solution[(length(w)+1) : (2*length(w))] ,1)
    } else{
        n2 <- rep(99999,length(w))
    }

    return(c(c2,n2))
  }

  # Search n1 if it is not given
  if(!is.numeric(n1)){
    #Score when using optimal second stage
    score <- function(n1){
      s2 <- stage_two(n1)
      n2 <- s2[ (length(s2)/2 + 1) : length(s2)]
      p <- score_direct(n1, n2)
      return(p)
    }

    optimum <- nloptr::nloptr(
      x0          = 82.33594,
      eval_f      = function(x) score(x),
      lb = 1,
      ub = 300,
      opts = list(
        # algorithm = "NLOPT_LN_NELDERMEAD",
        # NELDERMEAD failes sometimes. therefore, alternatively:
        algorithm = "NLOPT_LN_SBPLX",
        maxtime = 3600,
        xtol_rel = 0.01,
        #ftol_rel = 0.001,
        maxeval = 100
      )
    )
    n1 <- optimum$solution
  }


  s2 <- stage_two(n1)
  c2 <- s2[1 : (length(s2)/2)]
  v <- seq(-1, 3, 4/(length(c2)-1))

  n2_test <- s2[(length(s2)/2 + 1) : length(s2)]
  cf <-  max( 0.5 * (v[min(which(n2_test > 0))] + v[min(which(n2_test > 0)) - 1]) , -1)
  ce <-  min( 0.5 * (v[max(which(n2_test > 0))] + v[max(which(n2_test > 0)) + 1]) , 3)
  n2 <- n2_test[which(n2_test > 0)]

  dis <- (round(ce,1) - round(cf,1)) / (length(n2)-1)
  u <- seq(round(cf, 1), round(ce, 1), dis)

  n2_out <- function(z){
    spl <- splinefun(u, n2)
    p <- ifelse(cf <= z && z <= ce, spl(z) ,0)
    p <- p + lb_n2
    return(p)
  }

  c2_out <- function(z){
    spl <- splinefun(v, c2)
    p <- spl(z)
    return(p)
  }


  d <- design(cf, ce, n1, n2_out, c2_out)
  return(d)

}
