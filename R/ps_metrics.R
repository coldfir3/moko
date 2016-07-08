#' Distance betwen vector and matrix
#'
#' This function computes and returns the minimum distance between a vector and
#' a matrix
#'
#' #@export
#' @param point numeric vector
#' @param set numeric matrix
#' @param method the distance measure to be used. This must be one of
#'   "euclidean" or "manhattan" (default).
#' @return numeric value indicating the minimum distance between point and set
#' @examples
#' aps <- ps(matrix(rnorm(1:1000),ncol=2))
#' pdist(c(0,0), aps$set)
pdist <- function(point, set, method = "manhattan"){
  point <- t(as.matrix(point))
  dif <- t(apply(set, 1, function(set) set-point))
  if(method == "manhattan")
    d <- apply(dif, 1, function(dif) sum(abs(dif)))
  if(method == "euclidean")
    d <- apply(dif, 1, function(dif) sum((dif)^2)^0.5)
  return(min(d))
}

#' IGD: Inverted Generational Distance
#'
#' The IGD is a perfomance measure function of Pareto front fidelity and
#' corresponds to the average distance between all designs in the true set and
#' the closest design of the current set. Thus, the lower the IGD value, the
#' better the front is.
#'
#' \deqn{ \text{IGD}(\matx{T},\matx{P}) = \frac{1}{|\matx{T}|} \sum_{\vect{t}
#' \in \matx{T}} \text{min}(d(\vect{t} - \vect{p}))_{\vect{p} \in \matx{P}} }{}
#'
#' @references Shimoyama, K., Jeong, S., & Obayashi, S. (2013, June).
#'   Kriging-surrogate-based optimization considering expected hypervolume
#'   improvement in non-constrained many-objective test problems. In 2013 IEEE
#'   \emph{Congress on Evolutionary Computation} (pp. 658-665). IEEE.
#'
#' @param aps 'ps' object containing the 'actual' pareto front
#' @param tps 'ps' object containing the 'true' pareto front
#' @param norm Logical indicating if the fronts should be normalized (default =
#'   TRUE).
#' @inheritParams pdist
#'
#' @return returns the IGD metric
#' @export
#' @examples
#' aps <- ps(matrix(rnorm(1:1000),ncol=2))
#' tps <- ps(matrix(rnorm(1:2000),ncol=2))
#' igd(aps,tps)
igd <- function(aps, tps, method = "manhattan", norm = TRUE){
  if (norm){
    aps <- normalize(aps)
    tps <- normalize(tps)
  }
  n <- nrow(tps)
  d <- apply(tps, 1, function(point) moko:::pdist(point, aps, method))
  igd <- sum(d)/n
  return(igd)
}
