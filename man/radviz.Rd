% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualise.R
\name{radviz}
\alias{radviz}
\title{Plot a multiresponse or multivariate dataset indo a 2d radViz graph}
\usage{
radviz(data, ...)
}
\arguments{
\item{data}{data.frame containing the variables or observations to be ploted}

\item{...}{opitional plotting argumentos passed to \code{points} function.}
}
\description{
Description
}
\examples{
data <- data.frame(matrix(rnorm(1:50),ncol=5))
radviz(data, col='red')

}

