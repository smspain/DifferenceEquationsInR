## utility functions for analyzing difference equation models
## and other discrete dynamical systems.
## ----------------------------------------------------------------
## NOTES:
## Can I write a program to Taylor-linearize nonlinear models?
## Begin to package this with other specific linear map models and
## specific CA rules? Include Boolean networks?
##
## Last modified: 5 MAY 2016 (comments updated)
## Author: Seth M. Spain
### ---------------------------------------------------------------

## Functions for Boolean Arithmetic
#'@title Boolean not in modular arithmetic
#'@export
#'@examples
#' not(0)
#' not(1)
not <- function(x) (1 + x) %% 2

#' @title Addition modulo m. Defaults to Boolean arithmetic XOR.
#' @description A function to do addition modulo some integer, m.
#' @param x a number or a vector
#' @param y a number, defaults to NULL
#' @param m the modulus, defaults to 2 for binary addition (XOR)
#' @usage add.m(x,2) returns the sum mod 2 of the vector x
#' add.m(x,y,2) returns the (x + y) mod 2.
#' @export
add.m <- function(x,y=NULL,m=2) {
  if (length(x) > 1) {
    out <- sum(x) %% m
  } else {
    out <- (x + y) %% m
  }
  return(out)
}

#' @title Boolean arithmetic AND (multiplication mod 2)
#' @description Defaults to Boolean arithmetic AND, but can be
#' used for general multiplication mod m
#' @param x a number or a numeric vector
#' @param y a number, defaults to NULL
#' @param m the modulus, defaults to 2 for binary multiplication (AND)
#' @export
multiply.m <- function(x,y=NULL, m=2) {
    if (length(x) > 1) {
        out <- prod(x) %% m
    } else {
        out <- (x*y) %% m
    }
    return(out)
}

#' @title Boolean arithmetic OR. (addition mod 2)
#' @description (x + y + x*y) mod 2
#' @param x input from the set {0,1}/default operation
#' @param y input from the set {0,1}
#' @param m modulus for the operation
#' @export
#' @examples
#' or(0,0)
#' or(0,1)
#' or(1,0)
#' or(1,1)
or <- function(x,y, m=2) (x + y + x*y) %% m

#' @title Catalan numbers
#' @description calculate the nth Catalan number
#' @param n index of the Catalan number, e.g., n = 1, 1; n = 2, 2, n = 3, 5, ... etc.
#' @export
#' @examples catalan(n=4)
catalan <- function(n) choose(2*n,n)/(n+1)


#' @title Convert decimal integers to binary representation.
#' @description
#' Based on a function by Thomas Lumley
#' @param i an interger-valued number
#' @param type = ("vector", "binary"): should the results be returned
#' as a vector of binary digits or as a single binary number? Defaults
#' to vector

#' @examples
#' binary(90)
#' binary(150)
#' binary(150, type="binary")
#' @export
binary <- function(i, type="vector") {
  options(scipen=100, digits=4) # coerce to fixed notation, not scientific
  a <- 2^(0:9)
  b <- 2*a
  c <- sapply(i,function(x) sum(10^(0:9)[(x %% b)>=a]))
  if(type == "vector") {
      d <- strsplit(as.character(c),"")
      e <- as.numeric(unlist(d))
      return(e)
  } else {
     return(c)
  }
}


#' @title Converts a binary vector to its decimal equivalent
#' @description useful for calculating the Wolfram Code of a CA rule
#' @param bin A vector of length N, containing entries from the set {0,1}
#' with the least significant digit last (i.e., reads left to right)
#' @export
decimal <- function(bin) {
    N <- length(bin)
    x <- rev(bin)
    #cat(x,"\n")
    out <- sum(x*2^(0:(N-1)))
    return(out)
}

#' @title Hamming Distance calculator
#' @description Calculate the Hamming distance (# of locations where
#' two vectors differ)
#' @export
#' @examples
#' hamming.dist(c(0,0,0,0),c(0,0,1,0)) #should return 1
hamming.dist <- function(a,b) sum(ifelse(a == b, 1,0))

#' @title Calculate the "Up-Down" signature of a vector
#' @description The up-down signature of a vector (or sequence) is
#' the sign of the lag-1 forward differences of the vector
#' updown = sign(diff(vector))
#' @param x A numeric vector
#' @export
#' @examples
#' set.seed(12345)
#' x <- rnorm(10)
#' updown(x)
updown <- function(x) sign(diff(x))

#' @title "Graph paper" plot
#' @description generates a plot space with gridlines, like graph paper
#' @param xlim minimum and maximum values of the x-axis
#' @param ylim minimum and maximum values of the y-axis
#' @export
make.plot <- function(xlim=c(-5,5),ylim=c(-5,5)) {
  plot(NULL, xlim=xlim, ylim=ylim, axes=FALSE, type="n",
       xlab="", ylab="")
  for (x in xlim[1]:xlim[2] ) abline(v=x, col="gray", lty=3)
  for (y in ylim[1]:ylim[2] ) abline (h=y, col="gray", lty=3)
  axis(1, pos=0); axis(2, pos=0)
}

#' @title Find the middle value of a vector of size L
#' @description A helper function for initialization function
#' @param L length of the 1-D lattice
#' @export
middle <- function(L) {
    if (L %% 2 == 0) {
        #cat(L, "\n")
        out <- L/2
    } else {
        #(cat(L+1,"\n"))
        out <- (L+1)/2
    }
    return(out)
}

#' @title Initialize a vector for initial condition of a cellular automata simulation
#' @description Initializes a vector of a given size, with only one cell (the "middle") at 1,
#' all others 0; initial condition of the 1-D lattice.
#' @param L the length of the desired vector; corresponds to the size
#' of the one-dimensional lattice the cellular automata will "live" on
#' @param type which type of initialization, "middle" will set the center
#' cell to 1, all else to zero; "random" sets the lattice to a random configuration
#' with any cell equal to 1 with probability .5
#' @export
ca.initialize <- function(L, type=c("middle","random")) {
    if (type == "middle") {
     m <- middle(L)
     out <- c(rep(0,(m-1)),1,rep(0,(L-m)))
    } else if (type=="random") {
     out <- rbinom(n=L,size=1, prob = .5)
    } else {
        cat("Warning! type must be in c('middle', 'random')!", "\n")
    }
    return(out)
}



#' @title Prepare a space-time matrix
#' @description Prepare a matrix to store space-time evolution of
#' 1D cellular automata, Boolean network, or multiple variable dynamical
#' system
#' @param x0 a vector containing the initial state for all locations
#' @param nSteps number of iterations, defaults to 11
#' @export
## ----------------------------------------------------------------
prep.matrix <- function(x0, nSteps=11) {
  xout <- matrix(rep(NA, length(x0)*nSteps), nrow=nSteps)
  xout[1,] <- x0
  return(xout)
}

#' @title Rotate a Matrix for plotting
#' @description
#' adapted from cape package by
#' Anna L. Tyler and Wei Lu and Justin J. Hendrick and
#' Vivek M. Philip and Greg W. Carter
#' vectorized with help from SNAP Technologies at UAlaska-Fairbanks.
#' Primarily a utility function for plot.ca()
#' @keywords utility
#' @param a matrix (often containing the space-time diagram of a cellular
#' automaton model)
#' @export
rotate.mat <- function (mat) {
  new.mat <- t(mat)[,nrow(mat):1]
  rownames(new.mat) <- colnames(mat)
  colnames(new.mat) <- rownames(mat)
  return(new.mat)
}

#' @title Space-time plot
#' @description Generate the space time plot for 1D cellular automata models.
#' Useful for other applications, like percolation theory.
#' just a wrapper for the image function that calls the rotate.mat
#' function above and sets the margins neatly.
#' @param ca A matrix containing the space-time evolution for a 1D CA model
#' @param col Colours for the space-time plot (defaults to white for "off" and
#' black for "on")
#' @export
ca.plot <- function(ca, col=c("white", "black")) {
  opar <- par()
  par(mar=c(2,7,2,7)+0.1)
  image(rotate.mat(ca), axes=FALSE, col=col)
  suppressWarnings(par(mar=opar$mar))
}

#' @title calculate the greatest Liapunov exponent of a discrete dynamical system
#' @description Iterates the system forwards to a reasonable point, starting from
#' x0 and x0 + d0 (a small real number), then calculates the Liapunov exponent as
#' log(d.n/d.0) weighted by the inverse of the number of iterations between them.
#' ----------
#'
#' This is not a good operation to automate, use with extreme caution!!!
#'
#' @param x0 initial value
#' @param f the recursion relation describing the system (a function)
#' @param d0 the perturbation at time 0 (a small real number), defaults to .000001
#' @param nSteps how many iterations to move the system forward, defaults to 550
#' @param plot plot the two systems? defaults to FALSE
#' @return lambda the greatest Lyapunov exponent
#' @export
liapunov <- function(x0, f, d0=.000001, nSteps=550, plot=FALSE) {
    n <- nSteps
    k <- floor(n/2)
    x1 <- iterate.map(x0, f, nSteps=k)
    x0.d <- x1[k] + d0
    x1 <- iterate.map(x1[k], f, nSteps=(n-k))
    x2 <- iterate.map(x0.d,f,nSteps=(n-k))
    n <- length(x1)
    x1.n <- x1[n]
    xd.n <- x2[n]
    d.n <- xd.n - x1.n
    suppressWarnings(lambda <- (1/n)*log(abs(d.n)/abs(d0)))
    if (is.nan(lambda)) stop("Numerical error in log(d.n/d0)!")
    return(lambda)
}

#' @title Fit a linear model to a first-order difference equation
#' @description Function takes a vector of observed time-ordered data and fits a linear model for its
#' first order difference model
#' @param vec the vector of data
#' @param center set the first observation to zero and difference all observations from it, defaults to zero
#' @return the fit of the model, an object of class 'lm'
#' @export
lin.fit <- function(vec, center = FALSE) {
  x0 <- vec[1]
  N <- length(vec)
  if (center == TRUE) {
    x <- vec - x0
    before <- x[1:(N-1)]
    after <- x[2:N]
    fit <- lm(after ~ before)
  } else {
    before <- vec[1:(N-1)]
    after <- vec[2:N]
    fit <- lm(after ~ before)
  }
  return(fit)
}

#' @title Find the order of a polynomial equation generating a sequence of numbers
#' @description The function differences the observed sequence until the differenced sequence is
#' constant and returns the order of differencing and the value of the constant. For instance, if
#' the order of differencing returned is 2, then the sequence can be generated using a second order
#' polynomial. Constancy is determined by assessing whether the range of the difference
#' sequence is below a tolerance threshold (can be set to 0 for absolute constancy).
#' @param vec the input sequence
#' @param tol tolerance for determining whether an order of difference is "constant". Defaults to .01.
#' @param print logical, print the differenced sequences for each iteration, defaults to FALSE
#' @param diff.table logical, include a table of differences in output, defaults to FALSE
#' @examples
#' x <- c(1,3,6,10,15,21,28)
#' find.order(x)
#' @export
find.order <- function(vec, tol=.01, print=FALSE, diff.table=FALSE) {
  N <- length(vec)
  x <- vec
 if(diff.table == TRUE) {
  tab <- matrix(rep(NA, N*N), nrow=N)
 }
  order=0
  while(((range(x)[2] - range(x)[1]) <= tol) == FALSE) {
    x <- diff(x)
    if(print==TRUE) {cat("x = ", x, "\n")}
    if (diff.table==TRUE) {
      tab[(order+1),] <- fill(tab[(order+1),],x)
    }
    order=order + 1
  }
  if (diff.table==TRUE ) {
    tab <- rbind(vec, tab)
    }
  diff.val <- x[1]
  if (diff.table == TRUE) {
    out <- list("order"=order, "Difference value"=diff.val,"Difference Table"=tab)
  } else {
    out <- list("order"=order, "Difference value"=diff.val)
}
  return(out)
}

#' @title function to calculate the Nth triangular number
#' @description Returns the Nth triangular number (the sum of the natural numbers to N)
#' @param N the magnitude of the triangular number you want to calculate
#' @examples
#' naturals <- 1:15
#' for (n  in naturals ) {
#'   cat(triangular.n(n), "\n")
#'  }
#'  @export
triangular.n <- function(N) sum(1:N)

#' @title Runs test for randomness in a set of discrete values (especially binary)
#' @param y a sequence of (binary) data
#' @export
runs.test <- function(y) {
  n <- length(y)
  n.switch <- sum(y[2:n] != y[1:(n-1)])
  return(n.switch)
}

#' @title Return a sequence of Fibonacci numbers
#' @description Get the sequence of Fibonacci numbers to N. Defaults to the Fibonacci sequence starting
#' from 0: 0,1,1,2,3,5,8,13,21,...
#' Initial values default to 0 and 1, so the length of the sequence will be N+1. If the user provides
#' a starting value for y0 of 1, then the resulting sequence will be of length N. Both initial values
#' are user-settable, allowing exploration of any Fibonacci-type sequence. The function is a masked call to
#' the delta.two function, so users interested in more elaborate sequence control can use that to set the
#' parameters (the fibonacci function hard-sets them: a = 1, b = 1, c = 0).
#' @param N the last Fibonacci number in the sequence
#' @param y0 the first initial value, defaults to 0
#' @param y1 the second initial value, defaults to 1
#' @param monitor when TRUE prints the iterations starting from the second Fibonacci number. Defaults to
#' FALSE.
#' @examples
#' F <- fibonacci(N=15, monitor=TRUE)
#' F
#' @export
fibonacci <- function(N, y0=0,y1=1, monitor=FALSE) {
  if (y0==0) {
    N=N+1
  }
  out <- delta.two(y0=y0,y1=1,a=1,b=1,c=0,N=N, monitor=monitor)
  return(out)
}

#' @title Uses the Fibonacci sequence as the basis for triangle number summation.
#' @param N the final Fibonacci number to sum to
#' @export
fib.triangle <- function(N) {
  F <- fibonacci(N)
  y <- cumsum(F)
  return(y)
}

#' @title Find the real or complex roots of a quadratic equation
#' @description Find roots of the equation a*m^2 + b*m + c = 0
#' @param a the coefficient on the squared term
#' @param b the coefficient on the order-1 term
#' @param c the constant
#' @export
quad.roots <- function(a,b,c) {
    D = sqrt(as.complex(b^2 - 4*a*c))
    if (Im(D) == 0 ) {
        D = Re(D)
        root1 = (-b + D)/(2*a)
        root2 = (-b - D)/(2*a)
    } else {
        root1 = (-b + D)/(2*a)
        root2 = (-b - D)/(2*a)
    }
    return(list("Root.1"=root1,"Root.2"=root2))
}

#' @title Closed form solution for second order difference equation
#' initial value problem
#' @description Find solution to a second order difference equation:
#' s(k) = c1*m_1^k + c2*m_2^k, where m_1 and m_2 are the roots of the
#' characteristic equation for the second order difference equation.
#' @param root1 First root of the characteristic equation (from quad.roots)
#' @param root2 Second root of the characteristic equation
#' @param y0    time 0 initial value
#' @param y1    time 1 initial value
#' @param echo  defaults to TRUE, when TRUE,
#' prints  the closed form solution
#' @returns A list with the closed formula solution as a function of k,
#' constant c1, and constant c2
#' @export
closed.form2 <- function(root1,root2,y0,y1, echo=TRUE) {
    c1 = y0 - (y1 - root1)/root2
    c2 = y0 - c1
    solution = function(k) c1*root1^k + c2*root2^k
    if (echo == TRUE) {
        cat(c1,"*",root1,"^k + ",c2,"*",root2,"^k","\n" )
    }
    return(list("Solution"=solution,"c1"=c1, "c2"=c2))
}

