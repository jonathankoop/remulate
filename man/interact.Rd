% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effectsTie.R
\name{interact}
\alias{interact}
\title{interact}
\usage{
interact(param = NULL, indices, scaling = c("none", "std"))
}
\arguments{
\item{param}{numeric value or function with time parameter. Specifies the value of the effect for the statistic in the REM model}

\item{indices}{is a numeric vector of indices corresponding to the effects specified in \code{effects} argument of function \code{\link{remulateTie}} on which the interaction term needs to be computed.}

\item{scaling}{specifies the method for scaling the statistic after the interaction has been computed. \code{"none"} [default] gives raw value of the statistic at time t, \code{"std"} the statistic is standardized per time}
}
\description{
This function specifies the input for the send effect in the \code{formula} argument for the function \code{\link{remulateTie}} or \code{\link{remulateActor}}. Not to be used independently
}
