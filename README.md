# robustDAA
Robust Differential Abundance Analysis of Microbiome Sequencing Data

## Installation         

```
# install.packages(c("modeest", "quantreg", "lme4", "robustlmm", "lmerTest", "pbkrtest"))
# install.packages("devtools")
devtools::install_github("guanxunli/robustDAA")
```

### An Example
We illustrate the usage of robustDAA package using data included in the package.

```
## Load package and dataset
library(robustDAA)
## dta is a list consisting of Y as the OTU data set and Z as the metadata set.
Y <- dta$Y
dim(Y)
Z <- dta$Z
dim(Z)
formula <- "u"

## Run robustDAA  
res <- rlm_fun(
  otu_tab = Y, meta = Z, formula = paste("~", formula),
  reg_method = "psi.huber"
)
```
