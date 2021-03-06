% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StageTwo.R
\name{score_direct_smooth}
\alias{score_direct_smooth}
\title{Score version for smooth direct designs}
\usage{
score_direct_smooth(parameters, n1, n2, h, N, w)
}
\arguments{
\item{parameters}{Parameters specifying the design}

\item{n1}{First stage sample size}

\item{n2}{n_2-values}

\item{h}{Distance between two nodes}

\item{N}{4N+1 gives the number of nodes}

\item{w}{nodes inside the interval (cf,ce)}
}
\description{
\code{score_direct_smooth} gives the version of the score that is needed for \link{stage_two}.
}
