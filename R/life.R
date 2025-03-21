#' @title Create lattice for Conway's Game of Life
#' @description wrapper and extension of make.lattice. Appends cells to edges of the k x k lattice. The life
#' function assumes this sort of boundary (i.e., the finite lattice we simulate on isolated from the rest of
#' the theoretically infinite lattice by a stable boundary layer of cells that remain zero 'for all time').
#' @param k the size of the lattice (makes a square lattice), defaults to 10
#' @param prob the probability that a cell is occupied, defaults to .5
#' @export
make.life <- function(k=10, prob=.5) {
  L <- make.lattice(L=k, d=2)
  L <- site.perc(lattice=L, prob=prob)
  append <- rep(0,k)
  L <- rbind(append,L)
  L <- rbind(L,append)
  k=k+2
  append <- rep(0,k)
  L <- cbind(append,L)
  L <- cbind(L, append)
  return(L)
}


#' @title update function for Conway's Game of Life
#' @description Assumes a stable boundary that remains zero. Within this boundary, uses the Moore neighbourhood
#' around each cell to determine birth, life, and death:
#' (a) if a cell is currently occupied and two or three of its neighbouring cells are also occupied, the cell
#' remains alive, (b) if the cell is currently unoccupied and exactly three of its neighbours are occupied a
#' new "occupant" is born in the cell, and (c) if the cell is occupied and more than three of its 
#' neighouring cells are occupied, the occupant 'dies' - if the cell is unoccupied under this condition, it remains
#' empty.
#' @param L a lattice (typically provided by the make.life function)
#' @export
life <- function(L){
  k <- nrow(L)
  x <- L
  for (i in 2:(k-1)) {
    for (j in 2:(k-1)) {
      summa <- L[i-1,j-1] + L[i,j-1] + L[i+1,j-1] + L[i-1,j] + L[i+1,j] + L[i-1,j+1] + L[i,j+1] + L[i+1,j+1]
      if (x[i,j] == 1 & (summa == 2 | summa == 3)) {
        x[i,j] <- 1
      } else if (x[i,j] == 0 & summa == 3) {
        x[i,j] <- 1
      } else {
        x[i,j] <- 0
      }
    }
  }
  # for (i in 2:(k-1)) {
  #   summa1 <- L[i-1,k] + L[i,k] + L[1+1,k] + L[i-1,1] + L[i,1] + L[i+1,1] + L[i-1,2] + L[i,2] + L[i+1,2]
  #   summak <- L[i-1,k-1] + L[i,k-1] + L[1+1,k-1] + L[i-1,k] + L[i,k] + L[i+1,k] + L[i-1,1] + L[i,1] + L[i+1,1]
  #   if (x[i,1] == 1 & summa == 2 | summa == 3) {
  #     x[i,1] <- 1
  #   } else if (x[i,1] == 0 & summa == 3) {
  #     x[i,1] <- 1
  #   } else {
  #     x[i,1] <- 0
  #   }
  #   if (x[i,k] == 1 & summa == 2 | summa == 3) {
  #     x[i,k] <- 1
  #   } else if (x[i,k] == 0 & summa == 3) {
  #     x[i,k] <- 1
  #   } else {
  #     x[i,k] <- 0
  #   }
  # }
  # for (j in 2:(k-1)) {
  #   summa1 <- L[1,j]
  # }
  L <- x
  return(L)
}

#' @title iterate a 2-dimensional cellular automaton model
#' @description Takes a prepared lattice of cells and iterates a dynamical rule/update function to
#' produce a spacetime array for the lattice model. Outputs a three way array, such that the rows and columns
#' the cells of the lattice at a given time step such that each matrix slice contains the system state at
#' a given time point/iteration.
#' @param x the initial state matrix for the model
#' @param f the function defining the system's updates, the system's dynamical rule, defaults to 'life'
#' @param nSteps the number of time steps to iterate the system for
#' 
#' @export
iterate.2dca <- function(x, f=life, nSteps=10) {
  I <- ncol(x)
  J <- nrow(x)
  X <- array(rep(NA, I*J*nSteps), dim=c(I,J,nSteps))
  X[,,1] <- x
  for (n in 2:nSteps) {
    X[,,n] <- f(X[,,(n-1)])
  }
  return(X)
}