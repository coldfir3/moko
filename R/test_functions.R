#' Test function: The Nowacki Beam
#'
#' Description
#'
#' @param x vector of length 2 correspon the normalized beath and height of the beam
#' @param g vector of lenght 4 containing the limit of each constraint
#' @param l numeric length of the beam
#' @param F numeric force applied at the beam tip
#' @param E numeric elastic longitudinal moduli
#' @param G numeric elastic transversal moduli
#' @param v numeric poison ratio
#' @param box data.frame structure containing the upper and lower limits for \code{b} and \code{h}
#' @return vector of objective and constrain responses
#' @export
#' @examples
#'
#' nowacki_beam(c(0.5,0.5))
nowacki_beam <- function(x,
                         g = c(5, 240, 120, 10),
                         l = 1500, F = 5000,
                         E = 216620, G = 86650, v = 0.27,
                         box = data.frame(b = c(10, 50),h = c(20, 250))){

  b <- x[1] * (box$b[2] - box$b[1]) + box$b[1]
  h <- x[2] * (box$h[2] - box$h[1]) + box$h[1]

  Iy <- b*h^3/12
  Iz <- b^3*h/12
  J <- Iy + Iz

  A <- b*h
  S <- 6*F*l/(b*h^2)
  g1 <- F*l^3/(3*E*Iy) - g[1]
  g2 <- S - g[2]
  g3 <- 3*F/(2*b*h) - g[3]
  g4 <- h/b - g[4]

  return(c(A,S,g1,g2,g3,g4))
}

#' Test function: Shaffer1
#'
#' Description
#'
#' @export
#' @examples
#'
#' Shaffer1(c(0.5,0.5))
Shaffer1 <- function(x){
  f1 <- x^2
  f2 <- (x-2)^2
  return(c(f1,f2))
}

#' Test function: Shaffer2
#'
#' Description
#'
#' @export
#' @examples
#'
#' Shaffer2(c(0.5,0.5))
Shaffer2 <- function(x){
  f1 <- sapply(x,function(x){
    if(x<=1) return(x)
    else if (x<=3) return(x-2)
    else if (x<=4) return(4-x)
    else return(x-4)
  })
  f2 <- (x-5)^2
  return(c(f1,f2))
}

#' Test function: Fonseca and Fleming
#'
#' Description
#'
#' @export
#' @examples
#'
#' Fonseca(c(0.5,0.5))
Fonseca <- function(x){
  d <- length(x)
  f1 <- 1-exp(-sum(x-1/sqrt(d))^2)
  f2 <- 1-exp(-sum(x+1/sqrt(d))^2)
  return(c(f1,f2))
}

#' Test function: Kursawe
#'
#' Description
#'
#' @export
#' @examples
#'
#' Kursawe(c(0.5,0.5))
Kursawe <- function(x){
  d <- 3
  if (length(x) != d)
    stop('lenght(x) must be 3')
  f1 <- -10*(exp(-0.2*sqrt(x[1]^2+x[2]^2))+exp(-0.2*sqrt(x[2]^2+x[3]^2)))
  f2 <- abs(x[1])^0.8+5*sin(x[1]^3) + abs(x[2])^0.8+5*sin(x[2]^3) + abs(x[3])^0.8+5*sin(x[3]^3)
  return(c(f1,f2))
}

#' Test function: Viennet
#'
#' Description
#'
#' @export
#' @examples
#'
#' Viennet(c(0.5,0.5))
Viennet <- function(x){
  x <- x*6 - 3
  y <- x[2]
  x <- x[1]
  yy <- y^2
  xx <- x^2
  f1 <- 0.5*(xx+yy) + sin(xx+yy)
  f2 <- (3*x-2*y+4)^2/8 + (x-y+1)^2/27 + 15
  f3 <- 1/(xx+yy+1) - 1.1*exp(-(xx+yy))
  return(c(f1,f2,f3))
}

#' Test function: Viennet
#'
#' Description
#'
#' @export
#' @examples
#'
#' Binh(c(0.5,0.5))
Binh <- function(x){
  y <- x[2]
  x <- x[1]
  f1 <- 4*x^2 + 4*y^2
  f2 <- (x-5)^2 + (y-5)^2
  g1 <- (x-5)^2 + y^2 -25
  g2 <- - (x-8)^2 - (y+3)^2 + 7.7
  return(c(f1,f2,g1,g2))
}
