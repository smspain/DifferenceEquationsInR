## utility variables
.EMPTY    = 0
.OCCUPIED = 1
.FLOW     = 2

#' @title Make Lattice
#' @description generate a 1- or 2-dimensional square lattice
#' @export
make.lattice <- function(L,d=2) {
    if (d == 1) {
        out <- matrix(rep(NA,L), nrow=1)
    } else if (d == 2) {
        out <- matrix(rep(NA,L*L), nrow=L)
    } else {
        print("d must equal 1 or 2!")
    }
}

#' @title Proababilty a finite 1-d lattice will percolate
#' @description For an infinite lattice, the only site occupation probability that will
#' percolate on the lattice is 1. Any probability less than 1 will leave a proportion of
#' sites unoccupied that is equal to p.
#'
#' For a finite lattice, site occupation probabilities less than 1 can produce a lattice that percolates,
#' but the probability decreases with lattice length. For lattice length L and site occupation
#' probability p:
#'
#' Pr(percolates) = p^L
#' @param L an integer giving the length of the lattice or a vector of sites defining the lattice
#' @param p site occupation probability, defaults to .5
#' @examples perc.prob.1d(10)
#' perc.prob.1d(10, .55)
#' @export
perc.prob.1d <- function(L, p=.5) {
    if((length(L) > 1) == TRUE) {
        l <- length(L)
        perc.prob <- p^l
    } else  if(is.numeric(L) == TRUE) {
            perc.prob <- p^L
        } else {
        cat("Error! Input must be a natural number or a 1-dimensional lattice!", "\n")
    }
return(perc.prob)
}

#' @title Site percolation
#' @description fills the cells of a 1- or 2-dimensional square lattice
#' with probability given by 'prob'
#' @param lattice an initialized 1- or 2-dimensional lattice
#' @param prob probability that a cell is occupied
#' @export
site.perc <- function(lattice, prob) {
    dims <- length(dim(lattice))
    if (dims == 1 ) {
        L <- length(lattice)
        lattice <- rbinom(L,1,prob)
        return(lattice)
    } else {
        n <- ncol(lattice)
        L <- nrow(lattice)
        for (i in 1:n ) {
            for ( j in 1:L ) {
                lattice[i,j] <- rbinom(1,1,prob)
            }
        }
        return(lattice)
        }
}

flow <- function(g, i = NA, j = NA) {
    # -> Cycle through cells in top row
    #
    if (is.na(i) || is.na(j)) {
        for (j in 1:ncol(g)) {
            g = flow(g, 1, j)
        }
        return(g)
    }
    #
    # -> Check specific cell
    #
    if (i < 1 || i > nrow(g) || j < 1 || j > ncol(g)) return(g)
    #
    if (g[i,j] == .OCCUPIED || g[i,j] == .FLOW) return(g)

    g[i,j] = .FLOW

    g = flow(g, i+1, j)        # down
    g = flow(g, i-1, j)        # up
    g = flow(g, i, j+1)        # right
    g = flow(g, i, j-1)        # left

    g
}

percolates <- function(g) {
    g <- flow(g)
    for (j in 1:ncol(g)) {
        if (g[nrow(g), j] == FLOW) return(TRUE)
    }
    return(FALSE)
}
