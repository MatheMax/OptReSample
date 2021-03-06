% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_direct.R
\name{score_direct}
\alias{score_direct}
\title{Compute the score}
\usage{
score_direct(parameters, cf, ce, nodes, c2, n1, n2)
}
\arguments{
\item{parameters}{Parameters (alpha, power, standardized effect) whith which you want to build your design}

\item{cf}{Boundary for stopping for futility after the first stage}

\item{ce}{Boundary for stopping for efficacy after the first stage}

\item{nodes}{Nodes for the stage two functions}

\item{c2}{c_2-values that correspond to \code{nodes}}

\item{n1}{First stage sample size}

\item{n2}{n_2-values that correspond to \code{nodes}}
}
\description{
\code{score_direct} computes the objective score approximately if the \code{c_2}-function and the \code{n_2}-function
are approximated via a set of splines.
}
\details{
Note that nodes, c2 and n2 have to be of the same length.
The score is given as the expected sample size of the design under the alternative hypothesis.
The integrals are approximated via the trapezial rule.
}
