# polya <- function(size=60, init.prop=.5, nSteps=50) {
#     urn <- rep(NA,size)
#     urn <- rbinom(size,1,init.prop)
#     prop <- rep(NA, nSteps)
#     for ( s in 1:nSteps ) {
#         select <- sample(urn,1)
#         urn <- c(select,urn)
#         N <- length(urn)
#         y <- sum(urn)
#         prop[s] <- y/N
#     }
#     return(list("urn"=urn, "prop"=prop))
# }
#
# pout <- polya()
