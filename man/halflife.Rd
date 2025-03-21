% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time-to-value-functions.R
\name{halflife}
\alias{halflife}
\title{Find the half-life of a linear dynamical system}
\usage{
halflife(a)
}
\arguments{
\item{a}{the multiplicative constant of the system
(factor for geometric progression described by the system)}
}
\description{
Find the half-life of a linear dynamical system
}
\examples{
halflife(.5)
halflife(.75)
halflife(.99)
}
