fn1 <- function(time,gamma){
    require(Rcpp)
    .Call('fn1',c(time),c(gamma),PACKAGE='fnoise')
}
