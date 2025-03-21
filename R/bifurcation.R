#' @title Plot a bifurcation diagram
#' @description Currently only works for one parameter models. Plots the long-term
#' behavior of a dynamical system as a function of varying its parameter value.
#' @param parmin minimum parameter value, defaults to 1
#' @param parmax maximum parameter value, defaults to 4
#' @param f function defining the map, defaults to logistic equation
#' @param x0 initial value to begin trajectories from (constant over parameter
#'  variations)
#' @param niter number of iterations to assess "long-term" behavior of the map
#' @param length length of the vector of parameter values to be assessed
#' @param usepar if TRUE, use maximum of parmin or xmin for the lower limit of xlim
#' in the underlying call to plot()
#' @export
#' @examples
#' bifurcation() #produces the logistic map's bifurcation diagram
#' bifurcation(f=gauss4bifurc, parmin=-1, parmax=1, ylim=c(-1,1.5),
#' usepar=TRUE, xmin=-1, length=201, x0=.2)
bifurcation <- function(parmin=1, parmax=4, f=logistic.eq, x0=.30,niter=100,
                        length=301, usepar = FALSE,
                        xmin=0, ylim=c(-.1,1), cex=0.1, pch=19,
                        xlab="Parameter value", ylab=expression(x[n]))  {
    if (usepar == TRUE) {
        zed <- xmin
    } else {
        zed <- ifelse(parmin >= xmin, parmin, xmin)
    }
    plot(NULL, xlim = c(zed, parmax), ylim = ylim, xlab = xlab,
         ylab = ylab)
    a <- seq(parmin, parmax, length.out=length)
            # step through many values of the parameter

    for (z in 1:length(a)) {
        xl <- rep(NA,niter)
        xl[1] <- x0
        for (i in 2:niter) {
            xl[i] <- do.call(f,list(xl[i-1],a[z]))
        }
        uval <- unique(xl[40:niter])  # keep "typical" states after "a long run"
        points(rep(a[z], length(uval)), uval, cex = cex, pch = pch)
    }
}
