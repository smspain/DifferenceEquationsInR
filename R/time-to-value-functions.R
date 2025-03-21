###############################################################################
# Last updated 29 August 2019 (updated 'timeto' function)
# Author: Seth M. Spain
###############################################################################

#' @title Number of time periods to a specific value
#' @description How long does a linear dynamical system take to reach
#' a given value from a given initial condition?
#' @param x0 initial value
#' @param a multiplicative constant (factor of the geometric progression)
#' @param value value to be attained
#' @export
#' @examples time2value(100,.5,50) #1
#' time2value(200,.5,50) #2
time2value <- function(x0,a, value) round((log(value) - log(x0))/log(a))

#' @title Find the half-life of a linear dynamical system
#' @param a the multiplicative constant of the system
#' (factor for geometric progression described by the system)
#' @export
#' @examples halflife(.5)
#' halflife(.75)
#' halflife(.99)
halflife <- function(a) round(-log(2)/log(a))

#' @title Number of time periods for an affine dynamical system to be
#' within a given "error" of its equilibrium
#' @description Function calculates the equilibrium for an empirical or
#' theoretical affine model, and returns how many time periods it takes
#' to get within a small margin of the equilibrium (the argument "final.diff",
#' which defaults to .01).
#'
#' Uses the absolute value of difference between starting and equilibrium
#' value to allow for increasing or decreasing (or oscillatory) convergence
#' to the equilibrium. The "error" system, e[n+1] = a*e[n], where
#' e[n] = y[n] - y.star, will
#' always decrease towards zero in an affine system.
#'
#' Result returned should be rounded to nearest integer for interpretation!
#' @export
timeto <- function(model=NULL,a=NULL,b=NULL,y0,final.diff=.01) {
    if (is.null(model) == FALSE ) {
        a <- unname(coef(model)[2])
        b <- unname(coef(model)[1])
    } else {
        a <- a
        b <- b
    }
    ystar <- fixed.pt(a=a,b=b)
    out <- (log(final.diff) - log(abs(y0 - ystar)))/log(a)
    return(out)
}
