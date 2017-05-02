fn1 <- function(time,gamma){
    require(Rcpp)
    .Call('fn1',c(time),c(gamma),PACKAGE='fnoiser')
}

fn2 <- function(time,gamma,correlation){
    require(Rcpp)
    .Call('fn2',c(time),c(gamma),c(correlation),PACKAGE='fnoiser')
}
