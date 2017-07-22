#' AR1 Simulations and Estimations
#' 
#' @param n.draws   number of draws
#' @param n         length of the time series (in quarterts, incl 1 year foracast)
#' @param out_path  where to save the data
#' @examples
#' out_path <- system.file(package = "tdmc", "out")
#' # use n.draws = 1000 to replicate paper
#' sim_ar1(n.draws = 10, out_path = out_path)  
#' @export
#' @importFrom reshape2 melt
sim_ar1 <- function(n.draws = 10, n = 92, out_path = system.file(package = "tdmc", "out")){

  stopifnot(n.draws > 2)

  # --- Generate Data ------------------------------------------------------------

  # Beta Specification
  beta <- c(10, 0.5)

  # Error Term Specification
  u.m50  <- random_matrix(type = "arima.sim", n = n, n.draws = n.draws, 
                            model = list(ar = -0.50))
  u.p00  <- random_matrix(type = "arima.sim", n = n, n.draws = n.draws, 
                           model = list(ar = 0.0000000000001))
  u.p50  <- random_matrix(type = "arima.sim", n = n, n.draws = n.draws, 
                            model = list(ar = 0.50))
  u.p85  <- random_matrix(type = "arima.sim", n = n, n.draws = n.draws, 
                            model = list(ar = 0.85))
  u.x.p00  <- random_matrix(type = "arima.sim.integrated", n = n, n.draws = n.draws, 
                          model = list(ar = 0.0000000000001))


  u.x.p50  <- random_matrix(type = "arima.sim.integrated", n = n, n.draws = n.draws, 
                          model = list(ar = 0.50))
  u.x.p85  <- random_matrix(type = "arima.sim.integrated", n = n, n.draws = n.draws, 
                          model = list(ar = 0.85))

  # library(tsbox)
  # ts_plot(ts(u.m50 , start = 1990, frequency = 4))
  # ts_plot(ts(u.p00, start = 1990, frequency = 4))
  # ts_plot(ts(u.p85, start = 1990, frequency = 4))
  # ts_plot(ts(u.x.p00, start = 1990, frequency = 4))
  # ts_plot(ts(u.x.p50, start = 1990, frequency = 4))
  # ts_plot(ts(u.x.p85, start = 1990, frequency = 4))

  # Indicator Series Specification (Factor increases realative variance)
  X   <- 5 * random_matrix(type = "rw", n = n, n.draws = n.draws, 
                                  model = list(drift = 0.5))

  y.m50 <- true_value_matrix(X = X, U = u.m50, beta = beta)
  y.p00 <- true_value_matrix(X = X, U = u.p00, beta = beta)
  y.p50 <- true_value_matrix(X = X, U = u.p50, beta = beta)
  y.p85 <- true_value_matrix(X = X, U = u.p85, beta = beta)
  y.x.p00 <- true_value_matrix(X = X, U = u.x.p00, beta = beta)
  y.x.p50 <- true_value_matrix(X = X, U = u.x.p50, beta = beta)
  y.x.p85 <- true_value_matrix(X = X, U = u.x.p85, beta = beta)


  # --- Fixed Models -------------------------------------------------------------

  mod <- list()

  mod$clm50_m50 <- td_sim(Y = y.m50, X = X, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clm50_p00 <- td_sim(Y = y.p00, X = X, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clm50_p50 <- td_sim(Y = y.p50, X = X, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clm50_p85 <- td_sim(Y = y.p85, X = X, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clm50_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clm50_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clm50_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "chow-lin-fixed", fixed.rho = -0.50)

  mod$clp00_m50 <- td_sim(Y = y.m50, X = X, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp00_p00 <- td_sim(Y = y.p00, X = X, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp00_p50 <- td_sim(Y = y.p50, X = X, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp00_p85 <- td_sim(Y = y.p85, X = X, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp00_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp00_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp00_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "chow-lin-fixed", fixed.rho = 0.0)


  mod$clp50_m50 <- td_sim(Y = y.m50, X = X, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp50_p00 <- td_sim(Y = y.p00, X = X, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp50_p50 <- td_sim(Y = y.p50, X = X, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp50_p85 <- td_sim(Y = y.p85, X = X, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp50_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp50_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp50_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "chow-lin-fixed", fixed.rho = 0.50)


  mod$clp85_m50 <- td_sim(Y = y.m50, X = X, method = "chow-lin-fixed", fixed.rho = 0.85)
  mod$clp85_p00 <- td_sim(Y = y.p00, X = X, method = "chow-lin-fixed", fixed.rho = 0.85)
  mod$clp85_p50 <- td_sim(Y = y.p50, X = X, method = "chow-lin-fixed", fixed.rho = 0.85)
  mod$clp85_p85 <- td_sim(Y = y.p85, X = X, method = "chow-lin-fixed", fixed.rho = 0.85)
  mod$clp85_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "chow-lin-fixed", fixed.rho = 0.85)
  mod$clp85_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "chow-lin-fixed", fixed.rho = 0.85)
  mod$clp85_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "chow-lin-fixed", fixed.rho = 0.85)

  mod$lip00_m50 <- td_sim(Y = y.m50, X = X, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip00_p00 <- td_sim(Y = y.p00, X = X, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip00_p50 <- td_sim(Y = y.p50, X = X, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip00_p85 <- td_sim(Y = y.p85, X = X, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip00_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip00_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip00_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "litterman-fixed", fixed.rho = 0.0)

  mod$lip50_m50 <- td_sim(Y = y.m50, X = X, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip50_p00 <- td_sim(Y = y.p00, X = X, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip50_p50 <- td_sim(Y = y.p50, X = X, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip50_p85 <- td_sim(Y = y.p85, X = X, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip50_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip50_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip50_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "litterman-fixed", fixed.rho = 0.50)

  mod$lip85_m50 <- td_sim(Y = y.m50, X = X, method = "litterman-fixed", fixed.rho = 0.85)
  mod$lip85_p00 <- td_sim(Y = y.p00, X = X, method = "litterman-fixed", fixed.rho = 0.85)
  mod$lip85_p50 <- td_sim(Y = y.p50, X = X, method = "litterman-fixed", fixed.rho = 0.85)
  mod$lip85_p85 <- td_sim(Y = y.p85, X = X, method = "litterman-fixed", fixed.rho = 0.85)
  mod$lip85_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "litterman-fixed", fixed.rho = 0.85)
  mod$lip85_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "litterman-fixed", fixed.rho = 0.85)
  mod$lip85_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "litterman-fixed", fixed.rho = 0.85)



  mod$dym50_m50 <- td_sim(Y = y.m50, X = X, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dym50_p00 <- td_sim(Y = y.p00, X = X, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dym50_p50 <- td_sim(Y = y.p50, X = X, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dym50_p85 <- td_sim(Y = y.p85, X = X, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dym50_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dym50_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dym50_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "dynamic-fixed", fixed.rho = -0.50)

  mod$dyp00_m50 <- td_sim(Y = y.m50, X = X, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp00_p00 <- td_sim(Y = y.p00, X = X, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp00_p50 <- td_sim(Y = y.p50, X = X, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp00_p85 <- td_sim(Y = y.p85, X = X, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp00_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp00_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp00_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "dynamic-fixed", fixed.rho = 0.0)


  mod$dyp50_m50 <- td_sim(Y = y.m50, X = X, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp50_p00 <- td_sim(Y = y.p00, X = X, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp50_p50 <- td_sim(Y = y.p50, X = X, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp50_p85 <- td_sim(Y = y.p85, X = X, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp50_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp50_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp50_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "dynamic-fixed", fixed.rho = 0.50)


  mod$dyp85_m50 <- td_sim(Y = y.m50, X = X, method = "dynamic-fixed", fixed.rho = 0.85)
  mod$dyp85_p00 <- td_sim(Y = y.p00, X = X, method = "dynamic-fixed", fixed.rho = 0.85)
  mod$dyp85_p50 <- td_sim(Y = y.p50, X = X, method = "dynamic-fixed", fixed.rho = 0.85)
  mod$dyp85_p85 <- td_sim(Y = y.p85, X = X, method = "dynamic-fixed", fixed.rho = 0.85)
  mod$dyp85_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "dynamic-fixed", fixed.rho = 0.85)
  mod$dyp85_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "dynamic-fixed", fixed.rho = 0.85)
  mod$dyp85_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "dynamic-fixed", fixed.rho = 0.85)


  # --- Auto Models --------------------------------------------------------------

  mod$xcl_m50 <- td_sim(Y = y.m50, X = X, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xcl_p00 <- td_sim(Y = y.p00, X = X, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xcl_p50 <- td_sim(Y = y.p50, X = X, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xcl_p85 <- td_sim(Y = y.p85, X = X, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xcl_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xcl_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xcl_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "chow-lin-maxlog", truncated.rho = -1)

  mod$xclt_m50 <- td_sim(Y = y.m50, X = X, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xclt_p00 <- td_sim(Y = y.p00, X = X, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xclt_p50 <- td_sim(Y = y.p50, X = X, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xclt_p85 <- td_sim(Y = y.p85, X = X, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xclt_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xclt_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xclt_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "chow-lin-maxlog", truncated.rho = 0.2)

  mod$xcle_m50 <- td_sim(Y = y.m50, X = X, method = "chow-lin-minrss-ecotrim")
  mod$xcle_p00 <- td_sim(Y = y.p00, X = X, method = "chow-lin-minrss-ecotrim")
  mod$xcle_p50 <- td_sim(Y = y.p50, X = X, method = "chow-lin-minrss-ecotrim")
  mod$xcle_p85 <- td_sim(Y = y.p85, X = X, method = "chow-lin-minrss-ecotrim")
  mod$xcle_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "chow-lin-minrss-ecotrim")
  mod$xcle_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "chow-lin-minrss-ecotrim")
  mod$xcle_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "chow-lin-minrss-ecotrim")



  mod$xli_m50 <- td_sim(Y = y.m50, X = X, method = "litterman-maxlog", truncated.rho = -1)
  mod$xli_p00 <- td_sim(Y = y.p00, X = X, method = "litterman-maxlog", truncated.rho = -1)
  mod$xli_p50 <- td_sim(Y = y.p50, X = X, method = "litterman-maxlog", truncated.rho = -1)
  mod$xli_p85 <- td_sim(Y = y.p85, X = X, method = "litterman-maxlog", truncated.rho = -1)
  mod$xli_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "litterman-maxlog", truncated.rho = -1)
  mod$xli_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "litterman-maxlog", truncated.rho = -1)
  mod$xli_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "litterman-maxlog", truncated.rho = -1)

  mod$xlit_m50 <- td_sim(Y = y.m50, X = X, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlit_p00 <- td_sim(Y = y.p00, X = X, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlit_p50 <- td_sim(Y = y.p50, X = X, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlit_p85 <- td_sim(Y = y.p85, X = X, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlit_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlit_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlit_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "litterman-maxlog", truncated.rho = 0.2)

  mod$xlie_m50 <- td_sim(Y = y.m50, X = X, method = "litterman-minrss")
  mod$xlie_p00 <- td_sim(Y = y.p00, X = X, method = "litterman-minrss")
  mod$xlie_p50 <- td_sim(Y = y.p50, X = X, method = "litterman-minrss")
  mod$xlie_p85 <- td_sim(Y = y.p85, X = X, method = "litterman-minrss")
  mod$xlie_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "litterman-minrss")
  mod$xlie_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "litterman-minrss")
  mod$xlie_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "litterman-minrss")


  mod$xdy_m50 <- td_sim(Y = y.m50, X = X, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdy_p00 <- td_sim(Y = y.p00, X = X, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdy_p50 <- td_sim(Y = y.p50, X = X, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdy_p85 <- td_sim(Y = y.p85, X = X, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdy_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdy_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdy_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "dynamic-maxlog", truncated.rho = -1)

  mod$xdyt_m50 <- td_sim(Y = y.m50, X = X, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdyt_p00 <- td_sim(Y = y.p00, X = X, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdyt_p50 <- td_sim(Y = y.p50, X = X, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdyt_p85 <- td_sim(Y = y.p85, X = X, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdyt_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdyt_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdyt_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "dynamic-maxlog", truncated.rho = 0.2)

  mod$xdye_m50 <- td_sim(Y = y.m50, X = X, method = "dynamic-minrss")
  mod$xdye_p00 <- td_sim(Y = y.p00, X = X, method = "dynamic-minrss")
  mod$xdye_p50 <- td_sim(Y = y.p50, X = X, method = "dynamic-minrss")
  mod$xdye_p85 <- td_sim(Y = y.p85, X = X, method = "dynamic-minrss")
  mod$xdye_x.p00 <- td_sim(Y = y.x.p00, X = X, method = "dynamic-minrss")
  mod$xdye_x.p50 <- td_sim(Y = y.x.p50, X = X, method = "dynamic-minrss")
  mod$xdye_x.p85 <- td_sim(Y = y.x.p85, X = X, method = "dynamic-minrss")


  # --- Output -------------------------------------------------------------------

  ### raw Dataframe

  rdq <- data.frame(do.call(cbind, lapply(mod, function(x) x$errors.q)))
  rdy <- data.frame(do.call(cbind, lapply(mod, function(x) x$errors.y)))
  rdi <- data.frame(do.call(cbind, lapply(mod, function(x) x$insample)))

  # beta1
  rdb <- data.frame(do.call(cbind, lapply(mod, function(x) x$beta.hat[2, ] - beta[2])))

  # rho
  rdr <- data.frame(do.call(cbind, lapply(mod, function(x) x$rho.hat)))

  mdq <- melt(rdq, id = NULL)
  mdy <- melt(rdy, id = NULL)
  mdi <- melt(rdi, id = NULL)
  mdb <- melt(rdb, id = NULL)
  mdr <- melt(rdr, id = NULL)

  # split names
  splitnames <- function(x){
    nm <- matrix(unlist(strsplit(as.character(x$variable), split = "_")), 
                 ncol = 2, byrow=TRUE)
    data.frame(method = nm[, 1], process = nm[, 2], value = x$value)
  }

  dfq <- splitnames(mdq)
  dfy <- splitnames(mdy)
  dfi <- splitnames(mdi)
  dfb <- splitnames(mdb)
  dfr <- splitnames(mdr)

  save(dfq, dfy, dfi, dfb, dfr, n.draws, beta, file = file.path(out_path, "data_ar1_RW.RData"))
  
  TRUE
}

