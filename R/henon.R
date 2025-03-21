#' @title Henon map (a 2-D map displaying chaos)
#' @description Described in Henon, 1976. Default values come from this paper.
#' x = 1 - a*x^2 + y # "fold"
#' y = b*x           # "flip" and contract
#' @param x x-coordinate, defaults to 0
#' @param y y-coordinate, defaults to 0
#' @param a the a parameter (coefficient on x[n]^2 in the x equation), default is 1.4
#' @param b the b parameter (coefficient on x[n] in the y equation), default is 0.3
#' @param nSteps the number of iterations, default is 10000
#' @export
#' @examples
#' out <- henon()
#' henon.plot(out)
henon <- function(x=0,y=0,a=1.4,b=0.3, nSteps=10000) {
    for (n in 1:(nSteps-1)) {
        x[n+1] <- 1 - a*x[n]^2 + y[n]
        y[n+1] <- b*x[n]
    }
    out <- cbind("x_n"=x,"y_n"=y)
    return(out)
}

#' @title Plot the output of the Henon Map function
#' @param h a henon object
#' @export
#' @examples
#' out <- henon()
#' henon.plot(out)
henon.plot <- function(h) {
    N <- nrow(h)
    x <- h[,1]
    y <- h[,2]
    plot(x,y, cex=.25, main=paste("First", N, "iterates of the Henon map"))
}

