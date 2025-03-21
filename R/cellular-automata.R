#' @title Wolfram Rule 22
#' @description The cell is active in the next timestep if exactly one
#' of it and its left and right nearest neighbours is active.
#' @export
R22 <- function(sigma) {
    tmp <- sigma
    N <- length(sigma)
    for (i in 2:(N-1)) {
        if ((tmp[i-1] == 1 & tmp[i] == 0 & tmp[i+1] == 0 )|
            (tmp[i] == 1 & tmp[i-1] == 0 & tmp[i+1] == 0 ) |
            ( tmp[i+1]==1 & tmp[i-1] == 0 & tmp[i] == 0)) {
            sigma[i] <- 1
        } else {
            sigma[i] <- 0
        }
    }
    if ((tmp[N] == 1 & tmp[1] == 0 & tmp[2] == 0 )|
        (tmp[1] == 1 & tmp[N] == 0 & tmp[2] == 0 ) |
        ( tmp[2]==1 & tmp[1] == 0 & tmp[N] == 0)) {
        sigma[1] <-  1
    } else {
        sigma[1] <- 0
    }
    if ((tmp[N-1] == 1 & tmp[N] == 0 & tmp[1] == 0 )|
        (tmp[N] == 1 & tmp[N-1] == 0 & tmp[1] == 0 ) |
        ( tmp[1]==1 & tmp[N-1] == 0 & tmp[N] == 0)) {
        sigma[N] <-  1
    } else {
        sigma[N] <- 0
    }
    return(sigma)
}

#' @title Wolfram rule 30
#' @description for S(p,q,r) the next state function is
#'
#' p + (1 - 2p)(q + r - qr)
#' @export
R30 <- function(x) {
    L <- length(x)
    tmp <- x
    for (i in 2:(L-1)) {
        x[i] <- tmp[i-1] + (1 - 2*tmp[i-1])*(tmp[i] + tmp[i+1] - tmp[i]*tmp[i+1])
    }
    x[1] <- tmp[L] + (1 - 2*tmp[L])*(tmp[1] + tmp[2] - tmp[1]*tmp[2])
    x[L] <- tmp[L-1] + (1 - 2*tmp[L-1])*(tmp[L] + tmp[1] - tmp[L]*tmp[1])
    return(x)
}

#' @title Wolfram rule 90
#' @description Addition modulo 2 of the cell's left and right
#' nearest neighbours
#' @export
R90 <- function(x) {
    n <- length(x)
    tmp <- x
    for (i in 2:(n-1)){
        x[i] <- add.m(tmp[i-1],tmp[i+1])
    }
    x[1] <- add.m(tmp[n],tmp[2])
    x[n] <- add.m(tmp[n-1],tmp[1])
    return(x)
}

#' @title Wolfram rule 110
#' @description a "Class 4" cellular automata model
#'
#' For S(p,q,r), the update rule is q + r - qr - pqr
#' @export
R110 <- function(x) {
    tmp <- x
    L <- length(x)
    for (i in 2:(L-1)) {
        x[i] <- tmp[i] + tmp[i+1] - tmp[i]*tmp[i+1] - tmp[i-1]*tmp[i]*tmp[i+1]
    }
    x[1] <- tmp[1] + tmp[2] - tmp[1]*tmp[1] - tmp[L]*tmp[1]*tmp[2]
    x[L] <- tmp[L] + tmp[1] - tmp[L]*tmp[1] - tmp[L-1]*tmp[L]*tmp[1]
    return(x)
}

#' @title Wolfram rule 170
#' @description Rule 170
#' depends on initialize() and prep.matrix()
#' Acts as a left-shift register with periodic boundary condition
#' (i.e., the lattice is a one-dimensional ring)
#' @param x an initialized space-time cellular automata matrix
#' @export
R170 <- function(x) {
  L <- length(x)
  b <- x
  for (i in 1:(L-1)) {
    x[i] <- b[i+1]
  }
  x[L] <- b[1]
  return(x)
}

#' @title Wolfram rule 184
#' @description Long term behaviour of this rule varies qualitatively
#' as a function of its initial state (can display Class 1, 2, or 3 behaviour).
#'
#' S(p,q,r): qr + (1 - q)p
#' @export
R184 <- function(x) {
    L <- length(x)
    tmp <- x
    for (i in 2:(L-1)) {
        x[i] <- tmp[i]*tmp[i+1] + (1 - tmp[i])*tmp[i-1]
    }
    x[1] <- tmp[1]*tmp[2] + (1 - tmp[1])*tmp[L]
    x[L] <- tmp[L]*tmp[1] + (1 - tmp[L])*tmp[L-1]
    return(x)
}

#' @title Wolfram rule 254
#' @description Only a fully empty neighbourhood returns an inactive cell,
#' so almost all initializations land in the 'all-on' attractor. The only
#' initial condition that does not is the 'all-zero' condition, which results
#' in the repulsive "all-off" attractor. That is, any perturbation will
#' eventually lead to all cells being active.
#'
#' S(p,q,r): 1 - (1-p)(1-q)(1-r)
#' @export
R254 <- function(x) {
    L <- length(x)
    tmp <- x
    for (i in 2:(L-1)) {
        x[i] <- 1 - (1 - tmp[i-1])*(1 - tmp[i])*(1 - tmp[i+1])
    }
    x[1] <- 1 - (1 - tmp[L])*(1 - tmp[1])*(1 - tmp[2])
    x[L] <- 1 - (1 - tmp[L-1])*(1 - tmp[L])*(1 - tmp[1])
    return(x)
}

majority.rule <- function(x) {
    L <- length(x)
    tmp <- x
    for ( l in 2:(L-1) ) {
        x[l+1] <- sign(tmp[l-1] + tmp[l] + tmp[l+1])
    }
    x[1] <- sign(tmp[L] + tmp[1] + tmp[2])
    x[L] <- sign(tmp[L-1] + tmp[L] + tmp[1])
    return(x)
}

#' @title iterate a 1D Cellular automaton (or Boolean network or system of equations)
#' @description This function will take a dynamical rule and an initialized vector
#' and iterate it forward in time for nSteps iterations
#' @param x0 initial state of the system
#' @param f dynamical rule or "updating function". Defaults to Wolfram rule 170.
#' @param nSteps number of time steps to iterate the simulation forward;
#' defaults to 10
#' @examples x0 <- initialization(10)
#' x <- iterate.1dca(x)
#' @export
iterate.1dca <- function(x0, f=R170, nSteps=10) {
    x <- prep.matrix(x0, nSteps=nSteps)
    L <- length(x0)
    N <- nSteps
    for (j in 2:N) {
        x[j,] <- f(x[j-1,])
    }
    return (x)
}


# update.2D <- function(x,g, tau, v, neighborhood="vonNeumann") {
#     x <- cbind(rep(0,nrow(x)),x,rep(0, nrow(x)))
#     x <- rbind(rep(0,ncol(x)),x, rep(0,ncol(x)))
#     I <- nrow(x); J <- ncol(x)
#     tmp <- x
#     #cat("vector is ", x, "\n")
#     for ( i in 2:(I-1) ) {
#         for ( j in 2:(J-1) ) {
#             #cat("Iteration ", i, "\n")
#             if (neighborhood=="vonNeumann") {
#                 neighbors <- c(tmp[i-1,j],tmp[i+1,j], tmp[i,j-1], tmp[i, j+1])
#             }
#             #cat("neighbors are ", neighbors, "\n")
#             if (x[i,j] == 0 ) {
#                 n <- sum(neighbors[neighbors == 1])
#                 p.birth <- 1 - (1 - g)^n
#                 x[i,j] <- rbinom(1,1,p.birth)
#             } else {
#                 if (x[i,j] == 1 ) {
#                     m <- sum(ifelse(neighbors == 2,1,0))
#                     p.infected <- 1 - (1 - tau)^m
#                     x[i,j] <- ifelse(rbinom(1,1,p.infected)==1,2,1)
#                 } else {
#                     x[i,j] <- ifelse(rbinom(1,1,v)==1,0,2)
#                 }
#             }
#         }
#     }
#     tdat <- x[-1,]; tdat <- tdat[]
#     x <- x[-c(1,I),-c(1,J)]
#     return(x)
# }
