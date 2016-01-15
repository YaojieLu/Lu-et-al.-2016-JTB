# Optimal stomatal conductance under stochastic rainfall

This repository contains all the R code to run the analyses and generate the figures in this paper:
  
Yaojie Lu, Remko A. Duursma, Belinda E. Medlyn, **Optimal stomatal behaviour under stochastic rainfall**. Journal of Theoretical Biology, accepted.

# Instructions

To generate **all figures**, you should be able to do:

```r
source("run.R")
```

This generates figures in PDF in the subdirectory `output/figures`.

# Dependencies

The code will also attempt to install any missing packages. If you have problems with package dependencies, here is the list of packages you need to have (and their dependencies). Note that the "plotBy" package is only available from https://bitbucket.org/remkoduursma/plotby.

```
optimx,car,plotBy
```
