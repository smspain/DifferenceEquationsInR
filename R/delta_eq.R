#' @title Simulate from 1st order linear (affine) difference
#' Equations
#' @description
#' Difference equations of this form can correspond to either
#' arithmetic growth or geometric growth, depending on their
#' pattern of parameters (and initial values).
#' for A = 1, B = b: arithmetic growth
#' for A = a, B = 0: geometric growth (or, really, B = b)
#' @param y0 initial value for sequence, defaults to 1
#' @param a slope parameter for linear model, defaults to 2
#' @param b constant parameter for linear model, defualts to 0
#' @param N length of sequence, forward iterations = N-1, defaults to 5
#' @param monitor, print iteration number and next value? defaults to FALSE
#' @keywords difference equation
#' @export
#' @examples
#' x <- delta.eq()
#' series.plot(x)
delta.eq <- function(y0=1, a=2, b=0, N=5, monitor= FALSE) {
    y <- rep(NA, N)
    y[1] <- y0
    for (n in 1:(N-1) ) {
        y[n+1] <- a*y[n] + b
        if ( monitor == TRUE) {
            cat (n, y[n+1], "\n")
        }
    }
    return(y)
}

#' @title Function for simple affine map
#' @description Facilitates use of the iterate.map() function for
#' studying behavior of simple affine maps; alternative to delta.eq() function
#' @export
affine.map <- function(x, a=.75, b=1) a*x + b

#' @title Function for simple linear map
#' @description. Facilitates use of the iterate.map() function for
#' studying behavior of simple linear maps; alternative to delta.eq() function.
#' Wrapper for the affine.map() function with b set to 0.
#' @export
linear.map <- function(x, a=.5) affine.map(x, a=a, b=0)

#' @title Second-order linear autonomous difference equation
#' @description
#' defaults to Harrod-Domar-Hicks model of the growth of
#' national income (oscillatory convergence version)
#' as described in S. Goldberg (1958) "Introduction to difference
#' Equations" Chapter 0.
#' Takes 2 initial values, as needed to define a second-order recurrence
#' @param y0 first intial value, defaults to 2
#' @param y1 second initial value, defaults to 3
#' @param a coefficient on "next" time step of y (y[n+1]), defaults to 1
#' @param b coefficient on current time step of y (y[n]), defaults to -0.5
#' @param c constant, defaults to 1
#' @param N number of iterations (includes first two initial values)
#' @param monitor print iteration number and next step value of y, defaults to FALSE
#' @keywords second order difference equation
#' @export
#' @examples
#' y <- delta.two()
#' series.plot(y)
#' # The Fibonacci sequence is a 2nd-order recurrence, with these settings
#' F <- delta.two(y0=0,y1=1, a=1, b=1, c=0, N=11)
#' F
#' series.plot(F)
## --------------------------------------------------------------
delta.two <- function(y0=2,y1=3,a=1,b=-0.5,c=1,N=25, monitor= FALSE) {
    y <- rep(NA,N)
    y[1] <- y0
    y[2] <- y1
    for (n in 1:(N-2)) {
        y[n+2] <-  a*y[n+1] + b*y[n] + c
        if ( monitor == TRUE) {
            cat(n, y[n+2],"\n")
        }
    }
    return(y)
}
