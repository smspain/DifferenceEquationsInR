#' @title Numerical integration by Riemann Sum
#' @description Numerically integrate a function f using Riemann's method
#' of sums of rectangles. This gives the definite integral with lower
#' limit a and upper limit b.
#' @param f the function to integrate
#' @param a lower limit, defaults to 0
#' @param b upper limit of integration, defaults to 1
#' @param nSteps number of integration (delta-x) steps to take, defaults to 1000
#' @param type should the function be evaluated on the left or right
#' side, defaults to left
#' @export
riemann.sum <- function(f, a=0, b=1 , nSteps = 1000, type="left") {
    delta.x <- (b - a)/nSteps
    n <- 1:nSteps
    if (type=="left") {
        x <- c(a, a + n[-nSteps]*delta.x)
    } else {
        x <- a + n*delta.x
    }
    integral <- sum(delta.x*f(x))
    return(integral)
}

#' @title Basic numerical integration by Euler's method
#' @description Numerically solve a single OBE by numerical integration
#' using Euler's method:
#' dy/dx = f(x)
#' (y[n+1] - y[n])/dx = f(x) for small dx
#' y[n+1] = y[n] + dx*f(x)
#' @param y0 starting value for y defaults to 1
#' @param f the function to be integrated over
#' @param dx step size for numerical integration, defaults to .01
#' @param nSteps number of integration steps, defaults to 1000
#' @export
#' @examples
#' f <- function(x) 2.75*x*(1-x)
#' init <- runif(1)
#' sol <- euler(init, f=f)
#' plot(y ~ times, data=sol, type="l")
euler <- function(y0=1, f, dx=.01, nSteps=1000) {
    sx <- dx*nSteps
    n <- seq(0, sx, by=dx)
    y <- rep(NA, nSteps+1)
    y[1] <- y0
    for (k in 1:nSteps) {
        y[k+1] <- y[k] + dx*f(y[k])
    }
    return(list("y"=y,"times"=n))
}

#' @title Simulate a zombie apocalypse
#' @description based on Munz, Hudea, Imad, & Smith? When zombies attack!
#' Solve a system of ODEs, based on MATLAB code by Philip Munz (uses
#' Euler's method).
#' dx/dt = lim deltax/deltat as deltat -> 0
#' dx/dt = f(x) ~ deltax/deltat = f(x), for small enough deltat
#' deltax = deltat*f(x)
#' deltax = x[t+1] - x[t], so
#' x[t+1] = x[t] + dt*f(x):
#' Additionally plots output
#' adapted for R by Seth Spain, 15 MAR 2016
#' @param z0 - initial number of zombie infectees (at time 0), defaults to 1
#' @param alpha -  "zombie destruction" parm, defaults to .005
#' @param beta - "new zombie' parm, defaults to .0095
#' @param zeta - "zombie resurrection" rate, defaults to .0001
#' @param delta - background death rate, defaults to .0001
#' @param stime - stopping time, defaults to 5 (days)
#' @param dt - time step for numerical solution by Euler method, defaults to .01
#' @param N - human population at time 0, defaults to 500, considered 1000s
#' @export
#' @examples
#' zout <- zombies()
zombies <- function(z0=1, alpha=.005,beta=.0095,zeta=.0001, delta=.0001,
                    stime=5, dt=.01, N=500) {
    n <- stime/dt
    t <- rep(0, n+1)
    s <- rep(0, n)
    z <- rep(0, n)
    r <- rep(0, n)
    ###
    s[1] <- N
    z[1] <- z0
    r[1] <- 0
    timepts <- seq(0,stime, by=dt)
    ###
    i <- 1
    #cat(" ", "Susceptibles", "Zombies", "Removeds", "Sum", "\n")
    for (k in 1:(length(timepts)-1)) {
        s[i+1] <- s[i] + dt*(-beta*s[i]*z[i] - delta*s[i]) # birthrate = bkgd deathrate
        z[i+1] <- z[i] + dt*(beta*s[i]*z[i] - alpha*s[i]*z[i] + zeta*r[i])
        r[i+1] <- r[i] + dt*(alpha*s[i]*z[i] + delta*s[i] + zeta*r[i])
        sum <- s[i+1] + z[i+1] + r[i+1]
        #cat(i, s[i+1], z[i+1], r[i+1], sum, "\n")
        i <- i + 1
        if (s[i] < 0 | s[i] > N ) {
            print("Break S!")
            break
        }
        if (z[i] > N | z[i] < 0) {
            print("Break Z!")
            break
        }
        if (r[i] < 0 | r[i] > N) {
            print("Break R!")
            break
        }
    }
    plot(timepts[1:i], s[1:i], type="l", xlab="Time",
         ylab="Pop'n (1000s)", lwd=2)
    lines(timepts[1:i], z[1:i], lty=3, lwd=1.75)
    legend("topright", legend=c("Suscepties","Zombies"), lty=c(1,3),
           lwd=c(2,1.75), inset=.01)
    return(list("susceptibles"=s, "zombies"=z,"removed"=r,"times"=timepts))
}

#' @title Newton's method
#' @description Newton's method to iteratively find the roots of a function
#' For a function P(x), the algorithm is
#' x_[n+1] = x_n - (P(x_n)/P'(x_n)), where P'(x) is the first derivative
#' of P(x).
#' @param f the function to find a root
#' @param x0 the starting value to evaluate
#' @param iter number of iterations to carry out the algorithm
#' @export
#' @examples f <- function(x) x^2 - 2
#' newton(f=f, x0= 99, iter=14)

newton <- function(f, x0=2, iter=100) {
    r <- body(f)
    f.prime <- D(r, "x")
    tf <- function(y) y - f(y)/eval(f.prime, list(x=y))
    y <- iterate.map(a0 = x0, tf,nSteps = iter )
    #cat(y)
    yfinal <- y[iter]
    return(list("final"=yfinal, "iterations"=y))
}


#' @title Newton's method for extrapolating a sequence
#' @description The method repeatedly differences a given sequence until
#' the differencing returns a constant sequence. This indicates the order
#' of the model producing the sequence. By adding up the final
#' quasidiagonal, the next value of the sequence can be determined.
#' Returns a matrix containing original sequence and its differences,
#' the next value of the given sequence, and the order of the sequence.
#' Method is derived from Isaac Newton.
#' @param x the sequence to be extrapolated
#' @examples
#' x <-  c(1,2,4,8,16,31,57,99,163)
#' extrapolate(x)
#' @export
extrapolate = function(x) {
    y = x
    k = length(x)
    mymat = matrix(rep(NA, k*k), ncol=k)
    mymat[1,] = x
    i = 2
    while (length(unique(y)) > 1) {
        y = diff(y)
        p = length(y)
        pad = c(y,rep(NA,k-p))

        mymat[i,] = pad
        i = i+1
    }

    rank = i - 1
    m = rank
    n = p
    tally = 0
    if (rank == k) cat("Warning: Number of iterations equals sequence length.
  The recurrence might be a variant of 2^n. \n \n")


    while (m > 0 ) {
        while (n < k+1) {
            tally = tally + mymat[m,n]
            m = m - 1
            n = n + 1
            # cat("rank = ", rank, "n = ", n, "\n")
        }
    }
    rank = rank - 1 #number of differencings, not line of table
    return(list("Difference.Table"= mymat, "next.value" = tally, "order" = rank))
}
