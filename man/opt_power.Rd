% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power.R
\name{opt_power}
\alias{opt_power}
\title{Compute the power of a design}
\usage{
opt_power(d, parameters)
}
\arguments{
\item{d}{An object of class \code{design}.}

\item{parameters}{Parameters speficying the situation. See \link{parameters} for details.}
}
\description{
\code{opt_power} gives the power of a design on the desired point alternative.
}
\details{
It equals \link{real_power} when \code{effect = parameters$mualt}.
}
