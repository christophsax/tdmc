# tdmc

Reproducible Monte-Carlo Evaluation of Regression-Based Temporal Disaggregation
Methods.

This R package contains the data, the code and the content to build the paper:


To install code and data, including all required packages:

```r
# install.packages("devtools")   # if not installed
devtools::install_github("christophsax/tdmc")
```

To run the simulations (set `n.draws = 1000` to replicate the results in the
paper):

```r
library(tdmc)
out_path <- path.package("tdmc", "out")
sim_ar1(n.draws = 10, out_path = out_path) 
sim_sarima(n.draws = 10, out_path = out_path) 
```

To build the paper:

```r
setwd(out_path)

```

