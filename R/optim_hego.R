devtools::use_package("GPareto")
#' Title
#'
#' Description
#'
#' @param x a vector representing the input for which one wishes to calculate EHI
#' @param model object of class \code{mkm} containing the objectives and constrains
#' @param paretoFront optional object of class \code{ps} containing the actual Pareto set
#' @param minimization logical indicating if the EHVI is minimizing all objectives
#'  (\code{TRUE} by default) or maximizing all objectives (\code{FALSE})
#' @inheritParams GPareto::crit_EHI
#' @return The Constrained Expected Hypervolume Improvement at \code{x}.
#' @export
#' @examples
#' # ------------------------
#' # The Nowacki Beam
#' # ------------------------
#' n <- 50
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, nowacki_beam))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1:2, lower=rep(0.1,d)))
#' grid <- expand.grid(seq(0, 1, , 50),seq(0, 1, , 50))
#' ehvi <- apply(grid, 1, EHVI, model) # this computation may take some time
#' contour(matrix(ehvi, 50))
#' points(model@design, col=ifelse(model@feasible,'blue','red'))
#' points(grid[which.max(ehvi),], col='green', pch=19)
EHVI <- function(x, model, control = NULL){
  if (class(model) != 'mkm')
    stop('The class of "model" must be "mkm"')
  if (length(model@objective) == 1)
    stop('Incorrect Number of objectives. Must be more than 1')
  modelcontrol <- model@control
  if (is.null(control$minimization))
    control$minimization <- TRUE
  if (is.null(control$paretoFront))
    control$paretoFront <- ps(model@response[,model@objective], control$minimization)$set
  if (is.null(control$nb.samp))
    control$nb.samp = 50
  if (is.null(control$seed))
    control$seed = 42
  if (is.null(control$nb.samp))
    control$nb.samp = 50
  if (model@j == 0)
    probg <- 1
  else {
    model_g <- model@km[-model@objective]
    pred_g <- predict(list2mkm(model_g), data.frame(t(x)), modelcontrol)
    s_g <- pred_g$sd
    m_g <- pred_g$mean
    probg <- prod(pnorm(-m_g/s_g))
  }
  if (is.null(control$refPoint)){
    if (control$minimization)
      control$refPoint <- as.matrix(apply(control$paretoFront, 2, min))
    else
      control$refPoint <- as.matrix(apply(control$paretoFront, 2, max))
  }
  ehvi <- GPareto::crit_EHI(x, model@km[model@objective], control$paretoFront, control, modelcontrol$type)
 # ehvi <- 1
  return(ehvi*probg)
}

devtools::use_package("GenSA")
#' Title
#'
#' Description
#'
#' @inheritParams GenSA::GenSA
#' @inheritParams EHVI
#' @return Vector. The best set of parameters found.
#' @export
#' @examples
#' # ------------------------
#' # The Nowacki Beam
#' # ------------------------
#' n <- 20
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, nowacki_beam))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1:2, lower=c(0.1,0.1)))
#' max_EHVI(model)
max_EHVI <- function(model, lower = rep(0, model@d), upper = rep(1, model@d),
                     control = NULL, optimcontrol = list(max.time = 2)){
  if (class(model) != 'mkm')
    stop('The class of "model" must be "mkm"')
  fn <- function(x)
    return(-EHVI(x, model, control))
  res <- GenSA::GenSA(NULL, fn, lower, upper, control = optimcontrol)
  res$value <- -res$value
  return(res[c('value','par')])
}

#' Title
#'
#' Description
#'
#' @inheritParams GenSA::GenSA
#' @inheritParams EHVI
#' @return updated \link{\code{mkm}} model
#' @export
#' @examples
#' # ----------------
#' # Fonseca and Flemming
#' # ----------------
#' n <- 20
#' d <- 2
#' m <- 2
#' A <- 1 #verificar pq 4 nao funfa
#' fun <- Fonseca
#' doe <- 2*A*replicate(d,sample(0:n,n))/n - A
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res)
#' model <- hEGO(model, fun, 20, lower = -rep(A,d), upper = rep(A,d), quiet = FALSE)
#' tpf <- mco::nsga2(fun, d, 2, lower.bounds = -rep(A,d), upper.bounds = rep(A,d))$value
#' plot(tpf)
#' points(ps(model@response)$set, col = 'blue', pch = 19)
#' # ----------------
#' # Shaffer1
#' # ----------------
#' n <- 10
#' d <- 1
#' A <- 10
#' fun <- Shaffer1
#' doe <- 2*A*replicate(d,sample(0:n,n))/n - A
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res)
#' model <- hEGO(model, fun, 20, lower = -A, upper = A, quiet = FALSE)
#' tpf <- mco::nsga2(fun, d, 2, lower.bounds = -A, upper.bounds = A)$value
#' plot(tpf)
#' points(ps(model@response)$set, col = 'blue', pch = 19)
#' # ----------------
#' # Shaffer2
#' # ----------------
#' n <- 10
#' d <- 1
#' fun <- Shaffer2
#' doe <- 15 * replicate(d,sample(0:n,n))/n - 5
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res)
#' model <- hEGO(model, fun, 20, lower = -5, upper = 10, quiet = FALSE)
#' tpf <- mco::nsga2(fun, d, 2, lower.bounds = -5, upper.bounds = 10)$value
#' plot(tpf)
#' points(ps(model@response)$set, col = 'blue', pch = 19)
#' # ----------------
#' # Viennet
#' # ----------------
#' n <- 20
#' d <- 2
#' fun <- Viennet
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res)
#' model <- hEGO(model, fun, 80, quiet = FALSE)
#' pairs(ps(model@response)$set)
#' rgl::plot3d(ps(model@response)$set)
#' # ----------------
#' # Binh
#' # ----------------
#' n <- 20
#' d <- 2
#' fun <- Binh
#' doe <- cbind(rep(5,n), rep(3,n)) * replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objectives = 1:2))
#' model <- hEGO(model, fun, 40, upper = c(5,3), quiet = FALSE)
#' fun <- function(x) Binh(x)[c(1,2)]
#' gfun <- function(x) -Binh(x)[-c(1,2)]
#' tpf <- mco::nsga2(fun, d, 2, lower.bounds = c(0,0), upper.bounds = c(5,3),
#'                    constraints = gfun, cdim = 2)$value
#' plot(tpf)
#' points(ps(model@response[which(model@feasible),model@objective])$set, col = 'blue', pch = 19)
#' # ----------------
#' # The Nowacki Beam
#' # ----------------
#' n <- 20
#' d <- 2
#' fun <- nowacki_beam
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1:2, lower = rep(0.1,d)))
#' model <- hEGO(model, fun, 80, quiet = FALSE, control = list(rho = 0.1))
#' fun <- function(x) nowacki_beam(x)[c(1,2)]
#' gfun <- function(x) -nowacki_beam(x)[-c(1,2)]
#' tpf <- mco::nsga2(fun, d, 2, lower.bounds = c(0,0), upper.bounds = c(1,1),
#'                    constraints = gfun, cdim = 4)$value
#' plot(tpf)
#' points(ps(model@response[which(model@feasible),model@objective])$set, col = 'blue', pch = 19)
hEGO <- function(model, fun, nsteps, lower = rep(0, model@d), upper = rep(1, model@d), quiet = TRUE,
                 control = NULL, optimcontrol = list(max.time = 2)){
  time <- proc.time()
  if (class(model) != 'mkm')
    stop('The class of "model" must be "mkm"')
  for(n in 1:nsteps){
    x_star <- max_EHVI(model, lower, upper, control)$par
    y_star <- fun(x_star)
    model <- mkm(rbind(model@design,x_star), rbind(model@response,y_star), model@control)
    if (!quiet){
      cat('Current iteration:', n, '(elapsed', (proc.time()-time)[3], 'seconds)\n')
      cat('Current design:', round(x_star,3),'\n')
      cat('Current response:', round(y_star[model@objective],3),
          ifelse(tail(model@feasible,1),'(feasible)','(unfeasible)'),'\n\n')
    }
  }
  return(model)
}
