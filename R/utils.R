

random_matrix <- function(type = "arima.sim", n = 120, n.draws = 1, 
                         model = list(ar = 0.5), bootstrap = TRUE, 
                         no.neg = FALSE){
  
  # special case: calls standard case and integrates the result
  if (type == "arima.sim.integrated") {
    z.stationary <- random_matrix(type = "arima.sim", n = n, n.draws = n.draws, 
                                 model = model)
    z <- apply(z.stationary, 2, cumsum)
  } else {  # standard case
    z <- matrix(NA, nrow = n, ncol = n.draws)
    for (i in 1:n.draws){
      if (type == "arima.sim") {
        z.raw <- arima.sim(n = n, model = model, n.start = 1)
        # shorten if order[2] > 0
        if(!is.null(model$order[2])){
          if (model$order[2] > 0){
            z[, i] <- z.raw[-(1:model$order[2])]
          } else {
            z[, i] <- z.raw
          }
        } else {
          z[, i] <- z.raw
        } 
      } else if (type == "simulate.Arima") {
        z[, i] <- simulate(object = model, nsim = n, bootstrap = bootstrap)
      } else if (type == "rw") {
        z[, i] <- 500 + cumsum(model$drift + rnorm(n = n))
      } else{
        stop("type not defined")
      }
      if (no.neg){
        if (any(z[, i] < 0)){
          while (any(z[, i] < 0)){
            message(paste("Increasing series by 10000 in draw", i))
            z[, i] <- z[, i] + 10000
          }
        }
      }
    }
  }
  z
}

#' @export
#' 
true_value_matrix <- function(X, U, beta, random.beta = FALSE){
  Z <- X
  for (i in 1:dim(X)[2]){
    if (random.beta){
      Z[, i] <- runif(1) + colSums(runif(1) * t(X[, i]))
    } else {
      Z[, i] <- beta[1] + colSums(beta[-1] * t(X[, i]))
    }
  }
  Z + U
}




#' @export
#' @import tempdisagg
#' @import stats
#' 
td_sim <- function(Y, X, fr = 4, method, output = FALSE, ...){
  require(tempdisagg)
  # browser()
  if (!is.ts(Y)) Y <- ts(Y, frequency = fr)
  if (!is.ts(X)) X <- ts(X, frequency = fr)
  
  ### Aggregieren
  Y.long <- ta(Y, to = 1)
  Y.short <- window(Y.long, end = end(Y.long)[1] - 1)
  
  ### Schätzen
  Y.hat    <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  rho.hat  <- matrix(NA, nrow = 1, ncol = ncol(Y))
  beta.hat <- matrix(NA, nrow = 2, ncol = ncol(Y))
  R2  <- matrix(NA, nrow = 1, ncol = ncol(Y))

  for (i in 1:ncol(Y.short)){
    pb <- txtProgressBar(min = 1, max = ncol(Y.short), style=3)
    
    td.obj <- td(Y.short[, i] ~ 1 + X[, i], method = method, ...)
    Y.hat[, i]    <- predict(td.obj)
    rho.hat[, i]  <- td.obj$rho
    beta.hat[, i] <- coef(td.obj)[1:2]  # fix for dynamic
    R2[, i] <- td.obj$adj.r.squared
    cnames <- names(coef(td.obj))
    if (output){
      print(summary(td.obj))
    }
  
    setTxtProgressBar(pb, i)
    

  }

  dimnames(Y.hat) <- list(time(Y), paste("draw", 1:ncol(Y), sep = "_"))
  dimnames(rho.hat) <- list("rho", paste("draw", 1:ncol(Y), sep = "_"))
  dimnames(beta.hat) <- list(cnames[1:2], paste("draw", 1:ncol(Y), sep = "_"))
  
  z <- list()
  z$Y.hat    <- ts(Y.hat, fr = fr)
  z$Y <- ts(Y, fr = fr)
  
  z$Y.short <- Y.short
  z$rho.hat  <- as.numeric(rho.hat)
  z$beta.hat <- beta.hat
  z$R2 <- R2
  
  Y.y <- ta(z$Y)
  Y.hat.y <- ta(z$Y.hat)
  
  z$errors.y <- c(window(((Y.y) - (Y.hat.y)), start = end(Y.y)))
  
  errors <- (z$Y) - (z$Y.hat)
  z$errors.q <- c(window(errors, 
                         start = c(end(Y.y)[1], 1), end = c(end(Y.y)[1], 4)))
  
  z$insample <- c(window(errors, 
                         start = c(end(Y.y)[1]-1, 1), end = c(end(Y.y)[1]-1, 4)))


  z$fr <- fr
  z
}


  

PC <- function (x) {
  stopifnot(inherits(x,"ts"))
  z <- 100 * (x/L(x, -1) - 1)
  attributes(z) <- attributes(x)
  z
}

#' @export
PCY <- function(x){
  stopifnot(inherits(x,"ts"))
  z <- 100*(x/L(x, -frequency(x)) - 1)
  attributes(z) <- attributes(x) 
  z
}

#' @export
DIFF <- function(x, k = -1){
  stopifnot(inherits(x,"ts"))
  z <- x - L(x, k = k)
  attributes(z) <- attributes(x) 
  z
}

#' @export
DIFFY <- function(x){
  DIFF(x, k = -4)
}

#' @export
L <- function(x, k = -1){
  stopifnot(inherits(x, "ts"))
  z <- window(stats::lag(x, k), start = start(x), end = end(x), frequency = frequency(x), 
              extend = TRUE)
  attributes(z) <- attributes(x) 
  z
}

