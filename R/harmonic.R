#' @title Calculate the first through Nth Harmonic Numbers
#' @param N which Harmonic number do you want to get up to?
#' @export
#' @examples
#' harmonic()
harmonic <- function(N = 40){
    y <- rep(NA, N)
    y[1] <- 0
    for ( n in 1:(N-1) ){
        y[n+1] <- y[n] + 1/n
    }
    return(y)
}

