% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/c2.R
\name{c2}
\alias{c2}
\title{Computing optimal c_2-values}
\usage{
c2(parameters, z, n1, lambda1, lambda2, cf, ce)
}
\description{
\code{c2} returns the optimal value of the c_2-function in the Lagrangian framework
}
\details{
@param parameters Parameters specifying the design
  @param z Z_1-value of the first stage
  @param n1 First stage sample size
  @param lambda1 Penalization parameter for type I error
  @param lambda2 Penalization parameter for type II error
  @param cf Boundary for stopping for futility
  @param ce Boundary for stopping for efficacy
  @export
}
