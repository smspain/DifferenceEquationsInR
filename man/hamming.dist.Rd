% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/discrete-dynamic-utilities.R
\name{hamming.dist}
\alias{hamming.dist}
\title{Hamming Distance calculator}
\usage{
hamming.dist(a, b)
}
\description{
Calculate the Hamming distance (# of locations where
two vectors differ)
}
\examples{
hamming.dist(c(0,0,0,0),c(0,0,1,0)) #should return 1
}
