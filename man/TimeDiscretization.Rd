% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TimeDiscretization.R
\name{TimeDiscretization}
\alias{TimeDiscretization}
\title{function that discretized a vector of time with a given delta.}
\usage{
TimeDiscretization(rdata, subject, outcomes, predictors = NULL, Time,
  Delta)
}
\arguments{
\item{rdata}{input data which time point has to be discretized}

\item{subject}{subject identifiant}

\item{outcomes}{a vector of outcomes names}

\item{predictors}{independent variables  to be be included in the modeling}

\item{Time}{colname indicating the time}

\item{Delta}{discretized time step}
}
\value{
a discretized  time vector
}
\description{
function that discretized a vector of time with a given delta.
}
