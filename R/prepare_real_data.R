#' SARIMA Simulations from Real World Data templates
#' 
#' @param n.draws   number of draws
#' @param n         length of the time series (in quarterts, incl 1 year foracast)
#' @param bootstrap should the errors in the simulated series be bootstraped?
#                   (instead of drawn from a normal)
#' @examples
#' real_data <- prepare_real_data(10)  # use 1000 to replicate paper
#' @export
#' @import dynlm forecast
prepare_real_data <- function(n.draws = 10, n = 92, bootstrap = TRUE){

  # n.draws = 10
  # n = 92
  # bootstrap = TRUE

  set.seed(1)

  message("--- Load Data -----------------------------------------------------")

  data(template)

  assign_two_up <- function(x, dta){
    assign(x, na.omit(dta[, x]), envir = sys.frame(-3))
    NULL
  }
  sapply(colnames(data_m), assign_two_up, dta = ta(data_m, to = "quarterly"))
  sapply(colnames(data_q), assign_two_up, dta = data_q)


  message("--- Generate Data: CONSTRVZ ---------------------------------------")

    sbv_bt_hoch <- sbv_bt_hoch / 1000
    vz_bau <- vz_bau / 1000
    m.constrvz <- dynlm(sbv_bt_hoch ~ vz_bau)
    # summary(m.constrvz)
    m.constrvz.d <- (dynlm(L(log(sbv_bt_hoch), 4) ~ L(log(vz_bau), 4)))

    # Error Term Specification
    true.resid <- resid(m.constrvz)

    beta.constrvz <- coef(m.constrvz)

    # model of the true residuals
    # plot(true.resid)
    # auto.arima(true.resid)
    resid.mod <- Arima(true.resid, order = c(0, 0, 2), seasonal = c(0, 1, 0))

    # artificial residuals with the same model as the true residuals, 
    # normal innovations.
    u.constrvz  <- random_matrix(type = "simulate.Arima", 
                                n = n, n.draws = n.draws, 
                                model = resid.mod, bootstrap = bootstrap)

    # series.mod <- auto.arima(vz_bau)
    series.mod <- Arima(na.omit(vz_bau), order = c(1, 1, 1), seasonal = c(1, 0, 1))

    # Indicator Series Specification 
    X.constrvz   <- random_matrix(type = "simulate.Arima", 
                                 n = n, n.draws = n.draws, model = series.mod, 
                                 bootstrap = bootstrap)

    X.constrvz <- X.constrvz + 200

    y.constrvz <- true_value_matrix(X = X.constrvz, U = u.constrvz, beta = beta.constrvz)


  message("--- Generate Data: PHARMA -----------------------------------------")

    data("swisspharma")
    exports.q <- exports.q / 1000
    exports.q <- window(exports.q, 
                        start = start(sales.q), 
                        end = end(sales.q))
    m.pharma <- lm(sales.q ~ exports.q)
    m.pharma.d <- (lm(DIFFY(log(sales.q)) ~ DIFFY(log(exports.q))))

    # Beta Specification
    beta.pharma        <- coef(m.pharma)


    # Error Term Specification
    true.resid <- ts(resid(m.pharma), 
                     start = start(sales.q), 
                     frequency = frequency(sales.q))

    # summary(ur.kpss(true.resid, type="mu", lags="short"))

    # model of the true residuals
    # auto.arima(true.resid)
    # plot(true.resid)
    resid.mod <- Arima(true.resid, order = c(2, 0, 3), seasonal = c(1, 0, 1), include.mean=F)



    # artificial residuals with the same model as the true residuals, 
    # normal innovations.
    u.pharma  <- random_matrix(type = "simulate.Arima", 
                              n = n, n.draws = n.draws, 
                              model = resid.mod, bootstrap = bootstrap)


    # series.mod <- auto.arima(exports.q)
    series.mod <- Arima(exports.q, order = c(0, 1, 2), seasonal = c(0, 1, 1))

    # Indicator Series Specification 
    X.pharma   <- random_matrix(type = "simulate.Arima", 
                               n = n, n.draws = n.draws, model = series.mod, 
                               bootstrap = bootstrap, no.neg = TRUE)

    y.pharma <- true_value_matrix(X = X.pharma, U = u.pharma, beta = beta.pharma)


  message("--- Generate Data: IMP --------------------------------------------")

    m.imp <- dynlm(m06_chemie_n ~ m_t1)
    # summary(m.imp)
    # summary(dynlm(PCY(m06_chemie_n) ~ PCY(m_t1)))
    m.imp.d <- (dynlm(DIFFY(log(m06_chemie_n)) ~ DIFFY(log(m_t1))))

    # Beta Specification
    beta.imp        <- coef(m.imp)    # from actual data


    # Error Term Specification
    true.resid <- resid(m.imp)


    # model of the true residuals
    # plot(true.resid)
    # auto.arima(true.resid)
    resid.mod <- Arima(true.resid, order = c(1, 0, 2), seasonal = c(1, 0, 1), include.mean=F)

    # artificial residuals with the same model as the true residuals, 
    # normal innovations.
    u.imp  <- random_matrix(type = "simulate.Arima", 
                           n = n, n.draws = n.draws, 
                           model = resid.mod, bootstrap = bootstrap)

    # series.mod <- auto.arima(na.omit(m_t1))
    series.mod <- Arima(na.omit(m_t1), order = c(0,1,0), seasonal = c(1,0,0))

    # Indicator Series Specification 
    X.imp   <- random_matrix(type = "simulate.Arima", 
                            n = n, n.draws = n.draws, model = series.mod, 
                            bootstrap = bootstrap)

    while (any(X.imp < 20000)){
      X.imp <- X.imp + 10000
    }

    y.imp <- true_value_matrix(X = X.imp, U = u.imp, beta = beta.imp)


  message("--- Generate Data: LIK --------------------------------------------")
      
    m.lik <- dynlm(lik_dl ~ lik)
    # summary(m.lik)
    m.lik.d <- (dynlm(DIFFY(log(lik_dl)) ~ DIFFY(log(lik))))

    # Beta Specification
    beta.lik        <- coef(m.lik)    # from actual data


    # Error Term Specification
    true.resid <- resid(m.lik)


    # model of the true residuals
    # plot(true.resid)
    # auto.arima(true.resid)
    resid.mod <- Arima(true.resid, order = c(2, 0, 1), seasonal = c(0, 1, 1), include.mean=F)

    # artificial residuals with the same model as the true residuals, 
    # normal innovations.
    u.lik  <- random_matrix(type = "simulate.Arima", 
                           n = n, n.draws = n.draws, 
                           model = resid.mod, bootstrap = bootstrap)

    # series.mod <- auto.arima(na.omit(lik))
    series.mod <- Arima(na.omit(lik), order = c(0, 1, 2), seasonal = c(0, 0, 0))

    # Indicator Series Specification 
    X.lik   <- random_matrix(type = "simulate.Arima", 
                            n = n, n.draws = n.draws, model = series.mod, 
                            bootstrap = bootstrap, no.neg = TRUE)

    y.lik <- true_value_matrix(X = X.lik, U = u.lik, beta = beta.lik)


  message("--- Generate Data: ENERGY -----------------------------------------")


    # sapply(colnames(data_m), assign_globally, dta = data_m)


    m.energy <- dynlm(lik_energ ~ lik_oil)
    # summary(m.energy)
    m.energy.d <-  (dynlm(DIFFY(log(lik_energ)) ~ DIFFY(log(lik_oil))))


    # Beta Specification
    beta.energy        <- coef(m.energy)    # from actual data


    # Error Term Specification
    true.resid <- resid(m.energy)


    # model of the true residuals
    # plot(true.resid)
    # auto.arima(true.resid)
    resid.mod <- Arima(true.resid, order = c(1, 0, 1), seasonal = c(2, 0, 0), 
                       include.mean=F)

    # artificial residuals with the same model as the true residuals, 
    # normal innovations.
    u.energy  <- random_matrix(type = "simulate.Arima", 
                              n = n, n.draws = n.draws, 
                              model = resid.mod, bootstrap = bootstrap)

    # series.mod <- auto.arima(na.omit(lik_oil))
    series.mod <- Arima(na.omit(lik_oil), order = c(0, 1, 3), seasonal = c(0, 0, 0), include.drift = TRUE)

    # Indicator Series Specification 
    X.energy   <- random_matrix(type = "simulate.Arima", 
                               n = n, n.draws = n.draws, model = series.mod, 
                               bootstrap = bootstrap, no.neg = TRUE)


    y.energy <- true_value_matrix(X = X.energy, U = u.energy, beta = beta.energy)


  message("--- Collect Output ------------------------------------------------")

  z <- list()
  z$u.constrvz <- u.constrvz
  z$X.constrvz <- X.constrvz
  z$y.constrvz <- y.constrvz
  z$m.constrvz <- m.constrvz
  z$m.constrvz.d <- m.constrvz.d

  z$u.pharma <- u.pharma
  z$X.pharma <- X.pharma
  z$y.pharma <- y.pharma
  z$m.pharma <- m.pharma
  z$m.pharma.d <- m.pharma.d

  z$u.imp <- u.imp
  z$X.imp <- X.imp
  z$y.imp <- y.imp
  z$m.imp <- m.imp
  z$m.imp.d <- m.imp.d

  z$u.lik <- u.lik
  z$X.lik <- X.lik
  z$y.lik <- y.lik
  z$m.lik <- m.lik
  z$m.lik.d <- m.lik.d

  z$u.energy <- u.energy
  z$X.energy <- X.energy
  z$y.energy <- y.energy
  z$m.energy <- m.energy
  z$m.energy.d <- m.energy.d

  z

}
