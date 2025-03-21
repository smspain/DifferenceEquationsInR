#' @title calculations for linearizing the logistic model
#' @description Private function for calculating the coefficients
#' of a first-order Taylor series approximation of the logistic
#' equation
#' @export
lin.coef <- function(mu) {
    x <- logistic.fixedpt(mu=mu)
    f <- expression(mu*x*(1-x))
    f.prime <- expression(mu - 2*mu*x)
    alin <- eval(f.prime) #; alin
    blin <- eval(f) - eval(f.prime)*x #;blin
    return(list("alin"=alin,"blin"=blin))
}

#' @title Taylor Series linearization of the logistic equation
#' @description Approximates the logistic model near its equilibrium
#' @keywords approximation
#' @param y a value of the variable to be iterated
#' @param mu the parameter of the logistic model
#' @export
#' @examples
#' lin.logistic(y=.62,mu=2.5)
#' # for comparison
#' logistic.eq(x=.62, mu=2.5)
lin.logistic <- function(y, mu) {
    tmp <- lin.coef(mu)
    a <- tmp$alin
    b <- tmp$blin
    y.next <- a*y + b
    return(y.next)
}
