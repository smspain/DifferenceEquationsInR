#' @title Taylorize an arbitrary one-dimesional map
#' @description Linearize a nonlinear one-dimensional map using
#' a Taylor series expansion
#' @export
taylorize <- function(f, fixed.pt) {
    if(is.function(f)) {
        q <- body(f)
    } else if (is.expression(f)) {
        q <- f
    } else {
        cat("Error: f must be a function or an expression!")
    }
    x.star <- fixed.pt
    #cat(x.star, "\n")
    f.prime <- D(q, "x")
    #cat(f.prime, "\n")
    alin <- eval(f.prime, list(x=x.star))
    #cat(alin, "\n")
    blin <- eval(q, list(x=x.star)) - eval(f.prime, list(x=x.star))*x.star
    out <- function(y) alin*y + blin
    return(out)
    }
