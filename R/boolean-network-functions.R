#' @import igraph
# library(igraph)

#' @title Small example Boolean network
#' @description From Veliz-Cuba (2010) a four node Boolean network
#'
#' Note: This is not a Kauffman network, as nodes 1 and 2 have two inputs
#' each, while nodes 3 and 4 have one input each.
#'
#' x_1 = x_2^x_3
#'
#' x_2 = x_1^-x_4
#'
#' x_3 = x_2
#'
#' x_4 = -x_3
#' @param x an initial state vector for the 4 node network
#' @examples boolean.update(c(0,1,1,0))
#' @export
boolean.update <- function(x) {
    if(length(x) == 4) {
        tmp <- x
        x[1] <- tmp[2]*tmp[3]
        x[2] <- tmp[1]*not(tmp[4])
        x[3] <- tmp[2]
        x[4] <- not(tmp[3])
        return(x)
    } else {
        n <- length(x)
        cat("Length of state vector must be 4! Supplied vector is length", n, "\n")
    }
}

#' @title next state finder
#' @description Find next state for any input on a Boolean system
#' @param order the number of nodes in the Boolean system, defaults to 4
#' @param system the update rules for each node of the system, defaults to
#' 'boolean.update'
#' @export
next.states <- function(order=4, system=boolean.update) {
    states <- every.state(order)
    size <- 2^order
    transitions <- matrix
    transitions <- matrix(rep(NA, order*size), ncol=order)
    for(j in 1:size) {
        transitions[j, ] <- system(states[j,])
    }
    names <- c(paste("x",1:order,"_1", sep=""),paste("x",1:order,"_2", sep=""))
    out <- cbind(states,transitions)
    colnames(out) <- names
    return(out)
}

#' @title Every state for a Boolean system of order K
#' @description This function produces a matrix whose rows
#' contain every possible state for a Boolean system with K nodes.
#' @examples every.state(2)
#' every.state(4)
#' @export
every.state <- function(order) {
    k <- 2^order
    n <- 0:(k-1)
    m <- matrix(rep(0, order*k), ncol=order)
    for (j in 1:k ) {
        m[j,] <- fill(m[j,],binary(n[j]))
    }
    return(m)
}

#' @title Fill a vector with leading zeros
#' @description Add leading zeroes ahead of content to a vector
#' @param x an initialized vector (can be set to all NA) of length > length(a)
#' @param a the content (typically a vector with length less than length(x))
#' @export
fill <- function(x, a) {
    d <- length(x) - length(a)
    out <- c(rep(0,d),a)
    return(out)
}

#' @title Make the adjacency matrix for state transitions for a boolean network
#' @description Uses the next.states function to determine the state transition for all possible states of
#' a Boolean network of a given order and a given dynamic/updating rule. Intended to provide input
#' for graphing the state transitions (using igraph)
#' @examples
#' adjm <- make.adjacency()
#' g1 <- graph_from_adjacency_matrix (x)
#' plot(g1)
#' @export
make.adjacency <- function(order=4, system=boolean.update) {
    start.nextstate <- order+1
    end.nextstate <- order*2
    state.transitions <- next.states(order=order, system=system)
    K <- nrow(state.transitions)
    init.states <- 1:K - 0
    final.states <- rep(NA,K)
    transitions <- state.transitions[, start.nextstate:end.nextstate]
    for (k in 1:K) {
        final.states[k] <- decimal(transitions[k,])+1
    }
    adj <- matrix(rep(0,K*K), nrow=K)


    for (i in 1:K) {
        adj[i, final.states[i]] <- 1
    }
    rownames(adj) <- paste("S",1:16, sep="")
    colnames(adj) <- paste("S",1:16, sep="")
    return(adj)
}

#' @title Mask for igraph package's graph_from_adjacency_matrix
#' @description takes the output of the the make.adjacency function to produce a graph of the state transitions
#' for a given update rule on a set of variables. Intended for Boolean networks but will work with any
#' update function on a discrete and finite state space (though the graph could become unmanageably large.)
#' @export
transition.graph <- function(a, ...) igraph::graph_from_adjacency_matrix(a)


#' @title mask for the igraph package's plot function
#' @description provides ease-of-use for those unfamiliar with the igraph packe to easily plot state transition
#' graphs for finite state dynamical systems
#' @export
states.plot <- function(g, ...) igraph::plot.igraph(g)
