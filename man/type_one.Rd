% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TypeOne.R
\name{type_one}
\alias{type_one}
\title{Compute the type I error}
\usage{
type_one(parameters, cf, ce, nodes, c2)
}
\arguments{
\item{parameters}{Parameters (alpha, power, standardized effect) whith which you want to build your design}

\item{cf}{Boundary for stopping for futility after the first stage}

\item{ce}{Boundary for stopping for efficacy after the first stage}

\item{nodes}{Nodes for the c_2-function}

\item{c2}{Corresponding c_2-values}
}
\description{
\code{type_one} computes the type I error approximately if the \code{c_2}-function is approximated via a set of splines.
}
\details{
Note that nodes and c_2 have to be of the same length.
}
