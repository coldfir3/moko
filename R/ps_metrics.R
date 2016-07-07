#' Distance betwen vector and matrix
#'
#' This function computes and returns the minimum distance between a vector and a matrix
#'
#' #@export
#' @param point numeric vector
#' @param set numeric matrix
#' @param method the distance measure to be used. This must be one of "euclidean" or "manhattan" (default)
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

#' IGD pareto quality metric
#'
#' This function implements the method described in shimoyama2013kriging
#'
#' @inheritParams pdist
#' @param apf 'ps' object containing the 'actual' pareto front
#' @param tpf 'ps' object containing the 'true' pareto front
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
