# Reproducible Monte-Carlo Evaluation of Regression-Based Temporal Disaggregation Methods

The `tdmc` R package contains all *data*, *code* and *content* to build [*A Monte-Carlo Evaluation of Regression-Based Temporal Disaggregation Methods*](https://github.com/christophsax/tdmc/blob/master/inst/tex/document.pdf). Requires R and LaTeX, which are both free and open source.

To install, including all required R packages:

```r
# install.packages("devtools")   # if not installed
devtools::install_github("christophsax/tdmc")
```

To run the simulations:

```r
library(tdmc)
# Overwrites the existing simulation output (re-install to undo)
out_path <- system.file(package = "tdmc", "out")
sim_ar1(n.draws = 10, out_path = out_path) 
sim_sarima(n.draws = 10, out_path = out_path) 
```

Set `n.draws = 1000` to replicate the results in the paper. This will take 
about 2 hours on a laptop.

To build the paper:

```r
setwd(out_path)
library(knitr)
knit2pdf(system.file(package = "tdmc", "tex", "document.Rnw"))
```

To view the paper:
```r
browseURL("document.pdf")
```
