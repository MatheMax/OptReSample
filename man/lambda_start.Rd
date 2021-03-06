% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindLambda.R
\name{lambda_start}
\alias{lambda_start}
\title{Starting values for lambda}
\usage{
lambda_start(parameters)
}
\arguments{
\item{parameters}{Parameters specifying the design.}
}
\description{
\code{lambda_start} gives suitable starting values to find the values of lambda that yield to desired
type I - and type II - error.
}
\details{
These values are used as starting values in \link{find_lambda}.
The choice of the starting values seems quite arbitrarily. It is based on own examples. Improvments are very welcome.
}
