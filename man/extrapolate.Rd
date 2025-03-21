% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculus.R
\name{extrapolate}
\alias{extrapolate}
\title{Newton's method for extrapolating a sequence}
\usage{
extrapolate(x)
}
\arguments{
\item{x}{the sequence to be extrapolated}
}
\description{
The method repeatedly differences a given sequence until
the differencing returns a constant sequence. This indicates the order
of the model producing the sequence. By adding up the final
quasidiagonal, the next value of the sequence can be determined.
Returns a matrix containing original sequence and its differences,
the next value of the given sequence, and the order of the sequence.
Method is derived from Isaac Newton.
}
\examples{
x <-  c(1,2,4,8,16,31,57,99,163)
extrapolate(x)
}
