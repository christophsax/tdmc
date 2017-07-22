#' AR1 Simulations and Estimations
#' 
#' @param n.draws   number of draws
#' @param n         length of the time series (in quarterts, incl 1 year foracast)
#' @param out_path  where to save the data
#' @param bootstrap should the errors in the simulated series be bootstraped?
#'                  (instead of drawn from a normal)
#' @examples
#' out_path = path.package("tdmc", "out")
#' # use n.draws = 1000 to replicate paper
#' sim_sarima(n.draws = 10, out_path = out_path)  
#' @export
#' @importFrom reshape2 melt
sim_sarima <- function(n.draws = 10, 
                       n = 92, 
                       bootstrap = TRUE, 
                       out_path = path.package("tdmc", "out")){


  rd <- prepare_real_data(n.draws = n.draws, n = n, bootstrap = bootstrap)

  # --- Estimate -----------------------------------------------------------------

  mod <- list()

  mod$clm50_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clp00_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp50_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp85_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "chow-lin-fixed", fixed.rho = 0.85)

  mod$dym50_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dyp00_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp50_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp85_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "dynamic-fixed", fixed.rho = 0.85)

  mod$lip00_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip50_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip85_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "litterman-fixed", fixed.rho = 0.85)

  mod$xcl_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xclt_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xcle_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "chow-lin-minrss-ecotrim")

  mod$xdy_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdyt_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdye_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "dynamic-minrss")

  mod$xli_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "litterman-maxlog", truncated.rho = -1)
  mod$xlit_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlie_constrvz <- td_sim(Y = rd$y.constrvz, X = rd$X.constrvz, method = "litterman-minrss")


  mod$clm50_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clp00_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp50_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp85_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "chow-lin-fixed", fixed.rho = 0.85)

  mod$dym50_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dyp00_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp50_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp85_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "dynamic-fixed", fixed.rho = 0.85)

  mod$lip00_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip50_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip85_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "litterman-fixed", fixed.rho = 0.85)

  mod$xcl_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xclt_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xcle_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "chow-lin-minrss-ecotrim")

  mod$xdy_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdyt_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdye_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "dynamic-minrss")

  mod$xli_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "litterman-maxlog", truncated.rho = -1)
  mod$xlit_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlie_pharma <- td_sim(Y = rd$y.pharma, X = rd$X.pharma, method = "litterman-minrss")


  mod$clm50_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clp00_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp50_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp85_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "chow-lin-fixed", fixed.rho = 0.85)

  mod$dym50_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dyp00_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp50_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp85_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "dynamic-fixed", fixed.rho = 0.85)

  mod$lip00_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip50_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip85_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "litterman-fixed", fixed.rho = 0.85)

  mod$xcl_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xclt_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xcle_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "chow-lin-minrss-ecotrim")

  mod$xdy_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdyt_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdye_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "dynamic-minrss")

  mod$xli_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "litterman-maxlog", truncated.rho = -1)
  mod$xlit_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlie_imp <- td_sim(Y = rd$y.imp, X = rd$X.imp, method = "litterman-minrss")


  mod$clm50_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clp00_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp50_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp85_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "chow-lin-fixed", fixed.rho = 0.85)

  mod$dym50_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dyp00_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp50_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp85_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "dynamic-fixed", fixed.rho = 0.85)

  mod$lip00_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip50_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip85_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "litterman-fixed", fixed.rho = 0.85)

  mod$xcl_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xclt_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xcle_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "chow-lin-minrss-ecotrim")

  mod$xdy_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdyt_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdye_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "dynamic-minrss")

  mod$xli_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "litterman-maxlog", truncated.rho = -1)
  mod$xlit_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlie_lik <- td_sim(Y = rd$y.lik, X = rd$X.lik, method = "litterman-minrss")



  mod$clm50_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "chow-lin-fixed", fixed.rho = -0.50)
  mod$clp00_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "chow-lin-fixed", fixed.rho = 0.0)
  mod$clp50_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "chow-lin-fixed", fixed.rho = 0.50)
  mod$clp85_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "chow-lin-fixed", fixed.rho = 0.85)

  mod$dym50_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "dynamic-fixed", fixed.rho = -0.50)
  mod$dyp00_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "dynamic-fixed", fixed.rho = 0.0)
  mod$dyp50_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "dynamic-fixed", fixed.rho = 0.50)
  mod$dyp85_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "dynamic-fixed", fixed.rho = 0.85)

  mod$lip00_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "litterman-fixed", fixed.rho = 0.0)
  mod$lip50_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "litterman-fixed", fixed.rho = 0.50)
  mod$lip85_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "litterman-fixed", fixed.rho = 0.85)

  mod$xcl_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "chow-lin-maxlog", truncated.rho = -1)
  mod$xclt_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "chow-lin-maxlog", truncated.rho = 0.2)
  mod$xcle_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "chow-lin-minrss-ecotrim")

  mod$xdy_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "dynamic-maxlog", truncated.rho = -1)
  mod$xdyt_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "dynamic-maxlog", truncated.rho = 0.2)
  mod$xdye_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "dynamic-minrss")

  mod$xli_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "litterman-maxlog", truncated.rho = -1)
  mod$xlit_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "litterman-maxlog", truncated.rho = 0.2)
  mod$xlie_energy <- td_sim(Y = rd$y.energy, X = rd$X.energy, method = "litterman-minrss")
 

  # --- Output -------------------------------------------------------------------

  ### raw Dataframe

  rdq <- data.frame(do.call(cbind, lapply(mod, function(x) x$errors.q)))
  rdy <- data.frame(do.call(cbind, lapply(mod, function(x) x$errors.y)))
  rdi <- data.frame(do.call(cbind, lapply(mod, function(x) x$insample)))

  # rho
  rdr <- data.frame(do.call(cbind, lapply(mod, function(x) x$rho.hat)))

  mdq <- melt(rdq, id = NULL)
  mdy <- melt(rdy, id = NULL)
  mdi <- melt(rdi, id = NULL)
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
  dfr <- splitnames(mdr)

  save(dfq, dfy, dfi, dfr, n.draws, beta, file = paste0(out_path, "data_SARIMA.RData"))
  
  TRUE
}