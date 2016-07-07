devtools::use_package("DiceKriging")

#' @export
setClass('mkm', representation(
  km = 'list',
  objective = 'numeric',
  design = 'data.frame',
  d = 'numeric',
  n = 'numeric',
  m = 'numeric',
  j = 'numeric',
#  sd = 'numeric',
  response = 'data.frame',
#  norm_response = 'data.frame',
  feasible = 'logical',
  control = 'list'
  )
)

#' Multiobjective Kriging model
#'
#' This function creates a multiobjective kriging model. It is based on the \code{\link{km}}
#'  function of the \code{\link{DiceKriging}} package and creates a structured list of
#'  \code{\link{km}} objects.
#'
#' @param quiet logical indicating if the \code{print} commands inside \code{km} should be ommited
#' @param modelcontrol list of parameters passed to \code{\link{km}}.
#' @inheritParams DiceKriging::km
#' @return S4 class object that contains a list of \code{\link{km}} objects and also
#'  some usefull design data
#' @export
#' @examples
#' # ------------------------
#' # The Nowacki Beam
#' # ------------------------
#' n <- 10
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, nowacki_beam))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1:2))
mkm <- function(design, response, modelcontrol = NULL){

  design <- data.frame(design)
  names(design) <- paste('x',1:ncol(design),sep='.')
  response <- data.frame(response)
  names(response) <- paste('y',1:ncol(response),sep='.')

  model <- new('mkm')

  if (is.null(modelcontrol$objective))
    modelcontrol$objective <- 1:ncol(response)
  if (is.null(modelcontrol$quiet))
    modelcontrol$quiet <- TRUE
  if (is.null(modelcontrol$formula))
    modelcontrol$formula <- ~1
  if (is.null(modelcontrol$covtype))
    modelcontrol$covtype <- "matern5_2"
  if (is.null(modelcontrol$nugget.estim))
    modelcontrol$nugget.estim <- FALSE
  if (is.null(modelcontrol$estim.method))
    modelcontrol$estim.method <- "MLE"
  if (is.null(modelcontrol$optim.method))
    modelcontrol$optim.method <- "BFGS"
  if (is.null(modelcontrol$multistart))
    modelcontrol$multistart <- 1
  if (is.null(modelcontrol$gr))
    modelcontrol$gr <- TRUE
  if (is.null(modelcontrol$iso))
    modelcontrol$iso <- FALSE
  if (is.null(modelcontrol$scaling))
    modelcontrol$scaling <- FALSE
  if (is.null(modelcontrol$type))
    modelcontrol$type <- 'UK'
  if (is.null(modelcontrol$se.compute))
    modelcontrol$se.compute <- TRUE
  if (is.null(modelcontrol$cov.compute))
    modelcontrol$cov.compute <- FALSE
  if (is.null(modelcontrol$light.return))
    modelcontrol$light.return <- TRUE
  if (is.null(modelcontrol$bias.correct))
    modelcontrol$bias.correct <- FALSE
  if (is.null(modelcontrol$checkNames))
    modelcontrol$checkNames <- FALSE

  if (modelcontrol$quiet)
    invisible(capture.output(model@km <- apply(response, 2, function(response)
      DiceKriging::km(modelcontrol$formula,
                      design,
                      response,
                      modelcontrol$covtype,
                      modelcontrol$coef.trend,
                      modelcontrol$coef.cov,
                      modelcontrol$coef.var,
                      modelcontrol$nugget,
                      modelcontrol$nugget.estim,
                      modelcontrol$noise.var,
                      modelcontrol$estim.method,
                      modelcontrol$penalty,
                      modelcontrol$optim.method,
                      modelcontrol$lower,
                      modelcontrol$upper,
                      modelcontrol$parinit,
                      modelcontrol$multistart,
                      modelcontrol$control,
                      modelcontrol$gr,
                      modelcontrol$iso,
                      modelcontrol$scaling,
                      modelcontrol$knots,
                      modelcontrol$kernel))))
  else
    model@km <- apply(response, 2, function(response)
      DiceKriging::km(modelcontrol$formula,
                      design,
                      response,
                      modelcontrol$covtype,
                      modelcontrol$coef.trend,
                      modelcontrol$coef.cov,
                      modelcontrol$coef.var,
                      modelcontrol$nugget,
                      modelcontrol$nugget.estim,
                      modelcontrol$noise.var,
                      modelcontrol$estim.method,
                      modelcontrol$penalty,
                      modelcontrol$optim.method,
                      modelcontrol$lower,
                      modelcontrol$upper,
                      modelcontrol$parinit,
                      modelcontrol$multistart,
                      modelcontrol$control,
                      modelcontrol$gr,
                      modelcontrol$iso,
                      modelcontrol$scaling,
                      modelcontrol$knots,
                      modelcontrol$kernel))
  model@objective <- modelcontrol$objective
  model@design <- design
  model@response <- response
  model@d <- ncol(model@design)
  model@n <- nrow(model@design)
  model@m <- length(model@objective)
  model@j <- length(model@km) - model@m
#  model@sd <- (do.call(c,lapply(model@km, function(model) model@covariance@sd2)))^0.5
#  model@norm_response <- normalize(response)
  if(model@j == 0)
    model@feasible <- rep(TRUE, model@n)
  else
    model@feasible <- as.logical(apply(as.matrix(model@response[,-model@objective] <= 0), 1, prod))
  model@control <- modelcontrol
  return(model)
}

#' @export
setMethod("show", signature(object = "mkm"), function(object)
  print(cbind(object@design, object@response,feasible = object@feasible))
)

#' Predictor for a multiobjective Kriging model
#'
#' This functions performs predictions for a given dataset into a collection
#' of Kriging models (mkm object)
#'
#' @inheritParams DiceKriging::predict
#' @param object an object of class mpre
#'
#' @aliases predict
#'
#' @export
#' @examples
#' # ------------------------
#' # The Nowacki Beam
#' # ------------------------
#' n <- 100
#' d <- 2
#' N <- 50
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, nowacki_beam))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1:2, lower = rep(0.01, d)))
#' newx <- expand.grid(replicate(d,seq(0,1,,N),FALSE))
#' pred <- predict(model, newx)
#' realv <- t(apply(newx, 1, nowacki_beam))
#' par(mfrow=c(2,3), mar=c(2,2,1,1))
#' for (i in 1:6){
#'   contour(matrix((realv[,i]),N), col='red', lty=2, labels='')
#'   contour(matrix((pred$mean[,i]),N), add=T)
#' }
predict.mkm <- function(object, newdata, modelcontrol = NULL){
  if (is.null(modelcontrol))
    modelcontrol <- object@control
  output.list <- lapply(object@km, function(object)
    predict(object, newdata,
            modelcontrol$type,
            modelcontrol$se.compute,
            modelcontrol$cov.compute,
            modelcontrol$light.return,
            modelcontrol$bias.correct,
            modelcontrol$checkNames))
  output <- NULL
  for (i in 1:length(output.list)){
    output$trend <- cbind(output$trend, output.list[[i]]$trend)
    output$mean <- cbind(output$mean, output.list[[i]]$mean)
    output$sd <- cbind(output$sd, output.list[[i]]$sd)
    output$lower95 <- cbind(output$lower95, output.list[[i]]$lower95)
    output$upper95 <- cbind(output$upper95, output.list[[i]]$upper95)
  }
  output$norm_sd <- apply(normalize(output$sd), 1, prod)
  class(output.list) <- 'mpre'
  return(output)
}

#' @export
setMethod("predict", c("mkm"), predict.mkm)

#' @export
list2mkm <- function(list_of_models){ #rever e completar os outros slots
  model <- new('mkm')
  model@km <- list_of_models
  model@design <- as.data.frame(list_of_models[[1]]@X)
  model@response <- as.data.frame(lapply(list_of_models, function(model) model@y))
  model@d <- ncol(model@design)
  model@n <- nrow(model@design)
  model@m <- length(model@objective)
  model@j <- length(model@km) - model@m
  return(model)
}
