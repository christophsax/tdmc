#' LaTeX Regression Overviews of Template Regressions
#' 
#' @examples
#' regression_overview()
#' @import texreg
#' @export
regression_overview <- function(){
  rd <- prepare_real_data(1)
  z <- list()
  z$level <- texreg(list(rd$m.pharma,
            rd$m.imp, 
            rd$m.energy, 
            rd$m.constrvz, 
            rd$m.lik), 
       custom.model.names = c("SAL", "IMP", "ENE", "CON", "SER"),
       custom.coef.names = c("(Intercept)", 
            "Exports of Chemicals and Pharma", 
            "Total Imports", 
            "Consumer Prices Oil", 
            "Construction Employment", 
            "Consumer Prices Total"),
       booktabs = TRUE,
       caption = "OLS Regression in Levels ($x_t$)",
       label = "tab:reglevel",
       scriptsize = TRUE,
       caption.above = TRUE
       )    

  z$pc <- texreg(list(rd$m.pharma.d,
            rd$m.imp.d, 
            rd$m.energy.d, 
            rd$m.constrvz.d, 
            rd$m.lik.d), 
       custom.model.names = c("SAL", "IMP", "ENE", "CON", "SER"),
       custom.coef.names = c("(Intercept)", 
            "Exports of Chemicals and Pharma", 
            "Total Imports", 
            "Consumer Prices Oil", 
            "Construction Employment", 
            "Consumer Prices Total"),
       booktabs = TRUE,
       caption = "OLS Regression in seasonal Log-Differences ($\\log(x_t)-\\log(x_{t-4})$)",
       label = "tab:regdiff",
       scriptsize = TRUE,
       caption.above = TRUE
       )
  z
}
