## utility functions for analyzing difference equation models
## ----------------------------------------------------------------
## NOTES:
## Last modified: 29 August 2019 (updated fixed.pt function)
## Author: Seth M. Spain
### ---------------------------------------------------------------
## Begin functions for studying discrete maps
## begin with the standard affine map/1st order linear difference eq.
## first order, autonomous linear equation
## default value is "discrete exponential" growth
## --------------------------------------------------------------


#' @title Calculate steady state/equilibrium of a linear dyamical system
#' @description
#' fixed.pt() calculates the fixed point of a
#' first order, autonomous linear difference equation model.
#' Input can be either an lm-class model (model) OR a & b parameters
#' of a theoretical model
#' @param model a fitted model of class 'lm', defaults to NULL
#' @param a slope parameter of a linear recursion relation
#' @param b constant parameter of a linear recursion relation
#' @keywords equilibrium, fixed point
#' @export
#' @examples
#' fixed.pt(a=2,b=1)
fixed.pt <- function(model=NULL, a=NULL, b=NULL) {
  if (is.null(model) == FALSE) {
    a <- unname(coef(model)[2])
    b <- unname(coef(model)[1])
    if (a != 1) {
      y_star <- b/(1-a)
      return("Equilibrium"= y_star)
    } else {
      print ("a = 1! Equilibrium not defined!")
    }
  } else {
    if (a != 1) {
      y_star <- b/(1-a)
      return("Equilibrium" = y_star)
    } else {
      print ("a = 1! Equilibrium not defined!")
    }
  }
}

#' @title Solution to a first-order recurrence/difference equation
#' @description calculates solution  to
#'
#' y(n + 1) = r*y(n) + b
#'
#' in the form
#'
#' y(k) = r^k*(y0 - y.star) + y.star
#'
#' where y.star is the equilibrium, b/(1 - r).
#' @param x0 the initial condition for the system
#' @param model a least squares fit for the system (an lm object)
#' @param r as above
#' @param b as above
#' @export
first.solve <- function(x0, model=NULL, r=NULL, b=NULL) {
    if (is.null(model) == FALSE ) {
        r <- coef(model)[2]
        y.star <- fixed.pt(model)
        outf <- function(k) (r^k)*(x0 - y.star) + y.star
        return(outf)
    } else {
        y.star <- fixed.pt(a=r,b=b)
        outf <- function(k) (r^k)*(x0 - y.star) + y.star
        return(outf)
    }
}

#' @title Fixed point of a second order difference equation
#' @description Calculate the fixed point/equilibrium/steady state for
#' a 2nd order model of the form y[n+2] = a*y[n+1] + b*y[n] + c as
#' c/(1 - a - b). Can calculate from either a fitted model from lm or
#' user-specified parameters from a theoretical model.
#' @param model a fitted model of type lm defaults to NULL
#' @param a coefficient on y[n+1]
#' @param b coefficient on y[n]
#' @param c constant
#' @export
fixpt.2 <- function(model=NULL, a,b,c) {
    if ( is.null(model) == FALSE ) {
        a <- coef(model)[3]
        b <- coef(model)[2]
        c <- coef(model)[1]
    }
    out <- c/(1-a-b)
    names(out) <- "Equilibrium"
    return(out)
}

#' @title Calculate chi-square of a forward simulated empirical model
#' @description Chi-squared test statistic calculator for fit of
#' "forward simulated" model against observed data
#' @param data observed sequence of real values
#' @param sim the forward simulation sequence
#' @keywords chisq
#' @export
sim.chisq <- function(data, sim ) {
  N <- length(data)
  df <- N-2
  chi2 <- sum((sim[-1] - data[-1])^2/sim[-1])
  pval <- pchisq(chi2, df, lower.tail=FALSE)
  return(list("Chi.Sq" = chi2, "df"=df, "P-value" = pval))
}

#' @title Iterate a given user defined 1-dimensional map
#' @description
#' Allows modular examination of 1-dimensional discrete dynamical systems
#' @param a0 initial value of the sequence to be iterated
#' @param f "next value" function (f, which defaults to a function called "a")
#' @param nSteps the number of steps to iterate (include a0, forward iterations
#' equal to nSteps-1)
#' @keywords update, iterate
#' @export
#' @examples
#' iterate.map()
iterate.map <- function(a0, f=a, nSteps=11, ...) {
    x <- rep(NA, nSteps)
    x[1] <- a0
    for (n in 2:nSteps) {
        x[n] <- f(x[n-1])
    }
    return(x)
}


#' @title Identity function
#' @description useful utility for plotting cobweb analyses
#' identifies fixed points when identity(x) = map(x)
#' @param x identity function gives f(x) = x
#' @keywords identity
#' @export
identity <- function(x) x

#' @title Prep a cobweb analysis (graphical iteration)
#' @param f user defined function (update function, map), defaults to logistic
#' equation with parameter a = 2
#' @param from mask for curve function, defaults to 0
#' @param xlab, defaults to expression(x[n])
#' @param ylab, defaults to expression(x[n+1])
#' @keywords cobweb
#' @export
#' @examples
#' prep.cobweb()
prep.cobweb <- function(f=logistic.eq(x,mu=2), from=0, to=1, lwd=2,
                        col="darkgray", ylim=NULL,
                        xlab=expression(x[n]),
                        ylab=expression(x[n+1])) {
    curve(f, from=from, to=to, ylim=ylim, lwd=lwd, xlab=xlab,
          ylab=ylab)
    curve(identity(x), add=TRUE, lty=3, lwd=lwd, col=col)
}

#' @title Cobweb analysis of a 1-dimensional dynamical system
#' @description
#' "cobweb" a model (step through iterations); also called
#' "graphical iteration" (e.g., Sternberg, 2010).
#' @param x  a vector of the values of the sequence given by the difference
#' @param type defaults to segments, can also choose "arrows"
#' @param start.y initial position on y-axis, defaults to 0
#' @keywords cobweb
#' @export
#' @examples
#' cobweb()
cobweb <- function(x, type="segments", start.y=0, start.x=0, col="blue", lwd=2) {
  nsteps <- length(x)
  arrows(x[1],start.y,x[1],x[2], col=col, lwd=lwd)
  if(type == "segments") {
    for( i in 2:(nsteps-1) ){
      segments(x[i-1],x[i],x[i],x[i], col=col, lwd=lwd)
      #Sys.sleep(.05)
      segments(x[i],x[i],x[i],x[i+1], col=col, lwd=lwd)
      #Sys.sleep(.05)
    }
  } else {
    for( i in 2:(nsteps-1) ){
      arrows(x[i-1],x[i],x[i],x[i], col=col, angle=10, cex=.15, lwd=lwd)
      segments(x[i],x[i],x[i],x[i+1], col=col, lwd=lwd)
    }
  }
}


#' @title Plot the time series of a one-dimensional map
#' @description
#' x axis starts at value given by "start", ends at start + (length(x)-1)
#' @param x sequence given by the difference equation
#' @param start what value to start x-axis at, defaults to 1
#' @keywords orbit, trajectory
#' @export
#' @examples
#' plot.series()
series.plot <- function(x, start=1, xlab="Time period n",
                        ylab=expression(x[n]),
                        type="b", pch=19, lwd=1, lty=1, col="black", main="") {
    N <- length(x)
    end <- start+(N-1)
    plot(start:end,x,type=type, pch=pch,xlab=xlab, ylab=ylab, main=main,
         lwd=lwd, lty=lty, col=col)
}

#' @title Plot a trajectory/orbit for a 2-D map
#' @param trajectory, a 2 column matrix of the trajectory of a 2D Dynamical system
#' @keywords orbit, trajectory
#' @export
#' @examples
#' plot.orbit2D()
## -------------------------------------------------------------------
orbit2D <- function(trajectory,lty=1,pch=19,lwd=1,cex=1) {
  x <- trajectory[,1]
  y <- trajectory[,2]
  N <- length(x)
  plot(NULL, xlim=c(min(x)-.1, max(x)+.1),ylim=c(min(y)-.1, max(y)+.1),
       xlab="x", ylab="y")
  abline(v=0, lty=3, col="darkgray")
  abline(h=0, lty=3, col="darkgray")
  points(x[1],y[1],pch=pch,cex=cex)
  for (i in 1:(N-1) ) {
    segments(x[i],y[i],x[i+1],y[i+1],lty=lty,lwd=lwd)
    points(x[i+1],y[i+1], pch=pch)
  }
}


#' @title Implements the logistic equation
#' @description For use with the iterate.map() function. Function provides the 'next
#'  time step' value in the sequence {x[n]}. That is, for x[n] -> x[n+1]
#' @param x variable to iterate on
#' @param mu parameter of the logistic map x[n+1] = a*x[n]*(1 - x[n]), defaults to 2.5
#' @param N number of iterations (includes initial value), defaults to 20
#' @return for x = x[n], returns x[n+1]
#' @keywords logistic.eq
#' @export
#' @examples
#' logistic.eq(.5, a=2.5)
logistic.eq <- function(x, mu=2.5) mu*x*(1-x)

#' @title Calculate the fixed point of the logistic map:
#' @param a parameter of the logistic map x[n+1] = mu*x[n]*(1 - x[n])
#' @keywords equilibrium
#' @export
#' @examples
#' logistic.fixedpt()
logistic.fixedpt <- function(mu) 1 - (1/mu)

#' @title Tent map: A coordinate transform of the logistic map
#' @description Qualitatively similar behavior to the logistic map,
#' but very visually clear/intuitive.
#' Use with iterate.map() function:
#' out <- iterate.map(a0=init.value, f=tent.eq, nSteps=num.of.iterations)
#' @param x last value of x, i.e. x[n]
#' @keywords tent.map
#' @export
#' @examples tent.eq(.3)
tent.eq <- function(x) 2 - 2*abs(x)

#' @title Bernoulli shift map
#' @description funny little map, demonstratively useful, otherwise...?
#' @param x last value of x, i.e. x[n]
#' @keywords bernoulli.shift
#' @export
#' @examples
#' bernoulli.shift(.25)
#' bernoulli.shift(.5)
#' bernoulli.shift(2)
bernoulli.shift <- function(x) (2*x) %% 1

#' @title The Gauss map
#' @description x = exp(-a*x^2) + b \\
#' The Gauss iterated map is sometimes called the "mouse" map
#' because its bifurcation diagram looks like a mouse
#' @param x variable to iterate
#' @param a the a parameter, defaults to 4.9
#' @param b the b parameter, defaults to -0.58
#' @export
gauss.map <- function(x, a=4.9, b=-0.58) exp(-a*x^2) + b

#' @title Special version of the Gauss map for bifurcation plotting
#' @description gauss4bifurc holds the a-parameter of the Gauss map
#' constant at 4.9, while allowing the b-parameter to vary to produce
#' the "mouse" map.
#' @export
#' @examples
#' bifurcation(f=gauss4bifurc, parmin=-1, parmax=1, ylim=c(-1,1.5),
#' usepar=TRUE, xmin=-1, length=201, x0=.2)
gauss4bifurc <- function(x, b) gauss.map(x, a=4.9, b)

#' @title Discrete time SIR epidemiological model of infectious disease
#' @description A basic Susceptibles-Infectives-Removeds epidemiological model
#' for infectious disease. This is a compartment model defined by three equations
#' which describe the transfer of individuals between compartments: Susceptibles are
#' capable of becoming infected and enter the infectives compartment when they do.
#' Infectives enter the "removed" compartment by recovering or dying (both
#' are captured in a single rate parameter in this implementation.)
#' @param states the value of the S,I, and R compartments: should be done as percentage of population;
#' so initial S should be 1, common starting points for I might be .01 or .1, and R should be 0.
#' @param beta the transmission rate of the disease
#' @param alpha the rate defining transition from the I to the R compartment.
#' @param N the total population size to scale the beta coefficient, defaults to 1 so as to not
#' do any such scaling.
#' @export

discrete.SIR <- function(states, beta, alpha, N=1) {
  out <- rep(NA,length(states))
  out[1] <- states[1] - beta*(states[1]*states[2])/N
  out[2] <- states[2] + beta*(states[1]*states[2])/N - alpha*states[2]
  out[3] <- states[3] + alpha*states[2]
  if (out[1] < 0) out[1] <- 0
  if (out[2] < 0) out[2] <- 0
  if (out[3] < 0) out[3] <- 0
  if (out[1] > N) out[1] <- N
  if (out[2] > N) out[2] <- N
  if (out[3] > N) out[3] <- N
  return(out)
}
