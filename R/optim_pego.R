#' Tchebycheff
#'
#' Description
#'
#' @examples
#' grid <- expand.grid(seq(0, 1, , 50),seq(0, 1, , 50))
#'  res <- t(apply(grid, 1, nowacki_beam))
#'  plot(nowacki_beam_tps$x, xlim=c(0,1), ylim=c(0,1))
#'  grid <- grid[which(as.logical(apply(res[,-(1:2)] < 0, 1, prod))),]
#'  res <- res[which(as.logical(apply(res[,-(1:2)] < 0, 1, prod))),1:2]
#' for (i in 1:1000){
#' sres <- moko:::Tchebycheff(res[,1:2], s=100, rho=0.1)
#' points(grid[which.min(sres),], col='green')
#' }
Tchebycheff <- function(y, s=100, rho=0.1){ #add lambda as parameter or someway to define the wheighting
  lambda <- sample(0:s,ncol(y))/s
  lambda <- lambda/sum(lambda)
  y <- t(normalize(y))
 # y <- t(mco::normalizeFront(y))
  return((1-rho)*apply(lambda*y,2,max) + rho*apply(lambda*y,2,sum))
}

devtools::use_package("DiceOptim")
#' Constrained expected improvement
#'
#' description
#'
#' @export
#' @examples
#' # --------------------
#' # Branin-Hoo function
#' # --------------------
#' n <- 20
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' fun_cost <- DiceKriging::branin
#' fun_cntr <- function(x) 0.2 - prod(x)
#' fun <- function(x) return(cbind(fun_cost(x),fun_cntr(x)))
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1, lower=c(0.1,0.1)))
#' grid <- expand.grid(seq(0,1,,50),seq(0,1,,50))
#' ei <- apply(grid, 1, EI, model) # this computation may take some time
#' contour(matrix(ei,50))
#' points(model@design, col=ifelse(model@feasible,'blue','red'))
#' # ------------------------
#' # The Nowacki Beam
#' # ------------------------
#' fun_cost <- DiceKriging::branin
#' fun_cntr <- function(x){
#'  g1 <- 0.9 - sum(x)
#'  g2 <- sum(x) - 1.1
#'  g3 <- - x[1] + 0.5
#'  g4 <- x[2] - 0.5
#'  return(c(g1,g2,g3,g4))
#' }
#' fun <- function(x) return(c(fun_cost(x),fun_cntr(x)))
#' fun <- nowacki_beam
#' n <- 50
#' d <- 2
#' plugin <- NULL
#' doe <- replicate(d,sample(0:n,n))/n
#' doe <- replicate(d,sample(5:n,n-5))/n
#' for (i in 1:80){
#' res <- t(apply(doe, 1, fun))
#' res <- cbind(moko:::Tchebycheff(res[,1:2]),res[,-(1:2)])
#' #res <- cbind(res[,2],res[,-(1:2)])
#' model <- mkm(doe, res, modelcontrol = list(objective = 1, lower=rep(0.1,d)))
#' grid <- expand.grid(seq(0.1, 1, , 25),seq(0.1, 1, , 25))
#' ei <- apply(grid, 1, EI, model)
#' #, control=list(plugin=NULL)) # this computation may take some time
#' contour(matrix(ei, 25))
#' points(model@design, col=ifelse(model@feasible,'blue','red'))
#' points(grid[which.max(ei),], col='green', pch=19)
#' x_star <- unlist(grid[which.max(ei),])
#' y_star <- fun(x_star)
#' doe <- rbind(doe, x_star)
#' row.names(doe) <- NULL
#' }
EI <- function(x, model, control = NULL){
  if (class(model) != 'mkm')
    stop('The class of "model" must be "mkm"')
  if (model@m > 1)
    stop('Incorrect Number of objectives. Must be equal to 1')
  if (is.null(control$minimization))
    control$minimization <- TRUE
  if (model@j == 0)
    probg <- 1
  else {
    model_g <- model@km[-model@objective]
    pred_g <- predict(list2mkm(model_g), data.frame(t(x)), model@control)
    s_g <- pred_g$sd
    m_g <- pred_g$mean
    probg <- prod(pnorm(-m_g/s_g))
  }
  if (is.null(control$plugin)){
    if (control$minimization)
      control$plugin <- min(model@response[which(model@feasible),model@objective])
    else
      control$plugin <- max(model@response[which(model@feasible),model@objective])
  }
  ei <- DiceOptim::EI(x, model=model@km[[model@objective]],
                      plugin = control$plugin,
                      type = model@control$type,
                      minimization = control$minimization,
                      envir = control$envir)
  return(ei*probg)
}

devtools::use_package("GenSA")
#' Title
#'
#' Description
#'
#' @inheritParams GenSA::GenSA
#' @inheritParams EI
#' @return Vector. The best set of parameters found.
#' @export
#' @examples
#' # --------------------
#' # Branin-Hoo function
#' # --------------------
#' n <- 20
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- apply(doe, 1, DiceKriging::branin)
#' model <- mkm(doe, res)
#' max_EI(model)
max_EI <- function(model, lower = rep(0,model@d), upper = rep(1,model@d),
                   control = NULL, optimcontrol = list(max.time = 2)){
  if (class(model) != 'mkm')
    stop('The class of "model" must be "mkm"')
  if (model@m > 1)
    stop('model must have a single objective')
  fn <- function(x)
    return(-EI(x, model, control))
  res <- GenSA::GenSA(NULL, fn, lower, upper, control=optimcontrol)
  res$value <- -res$value
  return(res[c('value','par')])
}

#' MEGO: Multi-Objective Efficient Global Optimization Algorithm
#'
#' Executes \code{nsteps} iterations of the MEGO method to an object of class
#' \link{\code{mkm}}. At each step, a weighted kriging model is re-estimated (including
#' covariance parameters re-estimation) based on the initial design points plus
#' the points visited during all previous iterations; then a new point is
#' obtained by maximizing the Expected Improvement criterion (EI).
#'
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
#' model <- MEGO(model, fun, 20, lower = -rep(A,d), upper = rep(A,d), quiet = FALSE)
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
#' model <- MEGO(model, fun, 20, lower = -A, upper = A, quiet = FALSE)
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
#' model <- MEGO(model, fun, 20, lower = -5, upper = 10, quiet = FALSE)
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
#' model <- MEGO(model, fun, 80, quiet = FALSE)
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
#' model <- MEGO(model, fun, 40, upper = c(5,3), quiet = FALSE)
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
#' fun <- function(x) nowacki_beam(x, box = data.frame(b = c(10, 50),h = c(50, 250)))
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1:2, lower = rep(0.05,d)))
#' model <- MEGO(model, fun, 80, quiet = FALSE)
#' plot(model@design, col=ifelse(model@feasible,'blue','red'))
#' fun <- function(x) nowacki_beam(x, box = data.frame(b = c(10, 50),h = c(50, 250)))[c(1,2)]
#' gfun <- function(x) -nowacki_beam(x, box = data.frame(b = c(10, 50),h = c(50, 250)))[-c(1,2)]
#' tpf <- mco::nsga2(fun, d, 2, lower.bounds = c(0,0), upper.bounds = c(1,1),
#'                   constraints = gfun, cdim = 4)$value
#' plot(tpf)
#' points(ps(model@response[which(model@feasible),model@objective])$set, col = 'blue', pch = 19)
#'
#' # ----------------
#' # The Nowacki Beam (1D DEMO)
#' # ----------------
#' n <- 20
#' d <- 2
#' fun <- function(x) nowacki_beam(x)[c(1,3,6)]
#' doe <- replicate(d,sample(0:n,n))/n
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1, lower = rep(0.1,d)))
#' model <- MEGO(model, fun, 20, quiet = FALSE, control = list(rho = 0.1))
#' plot(model@design, col=ifelse(model@feasible,'blue','red'))
#'
#' #### SOME single objective optimization
#'
#' # -----------------------------------
#' # Branin-Hoo function (unconstrained)
#' # -----------------------------------
#' n <- 20
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' fun <- DiceKriging::branin
#' res <- apply(doe, 1, fun)
#' model <- mkm(doe, res, modelcontrol = list(lower=rep(0.1,d)))
#' model <- MEGO(model, fun, 20, quiet = FALSE)
#' plot(model@design, col=ifelse(model@feasible,'blue','red'))
#' # ---------------------------------------
#' # Branin-Hoo function (simple constraint)
#' # ---------------------------------------
#' n <- 10
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' fun_cost <- DiceKriging::branin
#' fun_cntr <- function(x) 0.2 - prod(x)
#' fun <- function(x) return(c(fun_cost(x),fun_cntr(x)))
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1, lower=rep(0.1,d)))
#' model <- MEGO(model, fun, 10, quiet = FALSE)
#' plot(model@design, col=ifelse(model@feasible,'blue','red'))
#' # ---------------------------------------
#' # Branin-Hoo function (narrow constraint)
#' # ---------------------------------------
#' n <- 10
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' fun_cost <- DiceKriging::branin
#' fun_cntr <- function(x){
#'  g1 <- 0.9 - sum(x)
#'  g2 <- sum(x) - 1.1
#'  g3 <- - x[1] + 0.5
#'  g4 <- x[2] - 0.5
#'  return(c(g1,g2,g3,g4))
#' }
#' fun <- function(x) return(c(fun_cost(x),fun_cntr(x)))
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1, lower=rep(0.1,d)))
#' model <- MEGO(model, fun, 10, quiet = FALSE)
#' plot(model@design, col=ifelse(model@feasible,'blue','red'))
#' # ---------------------------------------------
#' # Branin-Hoo function (discontinuos constraint)
#' # ---------------------------------------------
#' n <- 20
#' d <- 2
#' doe <- replicate(d,sample(0:n,n))/n
#' Griewank <-  function(x) {
#'  ii <- c(1:length(x))
#'   sum <- sum(x^2/4000)
#'   prod <- prod(cos(x/sqrt(ii)))
#'   y <- sum - prod + 1
#'   return(y)
#' }
#' fun_cost <- DiceKriging::branin
#' fun_cntr <- function(x) 1.6 - Griewank(x*10-5)
#' fun <- function(x) return(c(fun_cost(x),fun_cntr(x)))
#' res <- t(apply(doe, 1, fun))
#' model <- mkm(doe, res, modelcontrol = list(objective = 1, lower=c(0.1,0.1)))
#' model <- MEGO(model, fun, 20, quiet = FALSE)
#' plot(model@design, col=ifelse(model@feasible,'blue','red'))
MEGO <- function(model, fun, nsteps, lower = rep(0, model@d), upper = rep(1, model@d), quiet = TRUE,
                 control = NULL, optimcontrol = list(max.time = 2)){
  time <- proc.time()
  if (class(model) != 'mkm')
    stop('The class of "model" must be "mkm"')
  s_modelcontrol <- model@control
  if (is.null(control$s))
    control$s <- 100
  if (is.null(control$rho))
    control$rho <- 0.05

  design <- model@design
  response <- model@response
  if (model@m > 1){
    s_response <- moko:::Tchebycheff(response[,model@objective],
                                   s=control$s, rho=control$rho)
    s_response <- cbind(s_response, model@response[,-model@objective])
    s_modelcontrol$objective <- 1
  }
  else
    s_response <- model@response
  s_model <- mkm(design, s_response, s_modelcontrol)
  for(n in 1:nsteps){
    x_star <- max_EI(s_model, lower, upper, control, optimcontrol)$par
    y_star <- fun(x_star)
    design <- rbind(design, x_star)
    response <- rbind(response, y_star)
    rownames(response) <- NULL
    if (model@m > 1){
      s_response <- moko:::Tchebycheff(response[,model@objective],
                                       s=control$s, rho=control$rho)
      s_response <- cbind(s_response, response[,-model@objective])
    }
    else
      s_response <- response
    s_model <- mkm(design, s_response, s_modelcontrol)
    if (!quiet){
      cat('Current iteration:', n, '(elapsed', (proc.time()-time)[3], 'seconds)\n')
      cat('Current design:', round(x_star,3), '\n')
      cat('Current response:', round(y_star[model@objective],3),
          ifelse(tail(s_model@feasible,1),'(feasible)','(unfeasible)'),'\n\n')
    }
  }
  model <- mkm(design, response, model@control)
  return(model)
}

if(F){
n <- 20
d <- 2
grid <- expand.grid(replicate(d,seq(0,1,,50),F))
fun <- nowacki_beam
doe <- replicate(d,sample(0:n,n))/n
res <- t(apply(doe, 1, fun))

sresponse <- cbind(moko:::Tchebycheff(res[,1:2]),res[,-(1:2)])
smodel <- mkm(doe, sresponse, modelcontrol = list(objectives=1))
pred <- predict(smodel,grid)
s_g <- pred$sd[,-1]
m_g <- pred$mean[,-1]
probg <- apply(pnorm(-m_g/s_g), 1, prod)
contour(matrix(probg,ncol=50))
points(model@design, pch=19, col='blue')
(x_star <- max_EI(smodel, rep(0,d), rep(1,d), , list(max.time = 2)))
(y_star <- fun(x_star$par))
points(t(x_star$par), pch=19, col='green')
doe <- rbind(doe,x_star$par)
res <- rbind(res,y_star)
rownames(res) <- NULL
}
