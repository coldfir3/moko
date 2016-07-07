#'
#' Converts 'mco' objects into 'ps' objects
#' @param mco a 'mco' object
#' @return a 'ps' object
#'
mco2ps <- function(mco){
  if(all(class(mco)!='mco'))
    stop('object must be of "mco" class')
  ps <- NULL
  ps$set <- mco$value
  ps$x <- mco$par
  ps$index <- NULL
  ps$n <- nrow(ps$set)
  ps$m <- ncol(ps$set)
  class(ps) <- 'ps'
  return(ps)
}

#' @export
normalize <- function(data){
  ran <- apply(data, 2, range)
  data <- data.frame(t(apply(data, 1, function(data) (data - ran[1,])/(ran[2,]-ran[1,]))))
  return(data)
} # arrumar problema de NA (divizao por 0 )
