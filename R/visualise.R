#' Plot a multiresponse or multivariate dataset indo a 2d radViz graph
#'
#' Description
#'
#' @param data data.frame containing the variables or observations to be ploted
#' @param ... opitional plotting argumentos passed to \code{points} function.
#'
#' @export
#' @examples
#' data <- data.frame(matrix(rnorm(1:50),ncol=5))
#' radviz(data, col='red')
#'
radviz <- function(data, ...){

  k <- normalize(data)
  m <- ncol(data)
  n <- nrow(data)

  A.ang <- 2*pi*seq(0,m-1)/m
  A <- cbind(cos(A.ang), sin(A.ang))
  x <- apply(k, 1, function(k) sum(k*A[,1]))/apply(k, 1, sum)
  y <- apply(k, 1, function(k) sum(k*A[,2]))/apply(k, 1, sum)

  .pardefault <- graphics::par(no.readonly = TRUE)

  graphics::par(pty = "s", mar = c(1,1,1,1), bty="n",xaxt="n",yaxt="n")
  graphics::plot(A, xlab='', ylab='', pch=19, xlim=1.1*c(-1,1), ylim=1.1*c(-1,1))
  graphics::polygon(A, lty=3)
  graphics::text(1.1*A, names(data))
  graphics::points(x,y, ...)

  graphics::par(.pardefault)
}
