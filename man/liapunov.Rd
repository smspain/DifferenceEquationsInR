% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discrete-dynamic-utilities.R
\name{liapunov}
\alias{liapunov}
\title{calculate the greatest Liapunov exponent of a discrete dynamical system}
\usage{
liapunov(x0, f, d0 = 1e-06, nSteps = 550, plot = FALSE)
}
\arguments{
\item{x0}{initial value}

\item{f}{the recursion relation describing the system (a function)}

\item{d0}{the perturbation at time 0 (a small real number), defaults to .000001}

\item{nSteps}{how many iterations to move the system forward, defaults to 550}

\item{plot}{plot the two systems? defaults to FALSE}
}
\value{
lambda the greatest Lyapunov exponent
}
\description{
Iterates the system forwards to a reasonable point, starting from
x0 and x0 + d0 (a small real number), then calculates the Liapunov exponent as
log(d.n/d.0) weighted by the inverse of the number of iterations between them.
----------

This is not a good operation to automate, use with extreme caution!!!
}
