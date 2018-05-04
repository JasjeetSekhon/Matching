## Matching: Multivariate and Propensity Score Matching Software for Causal Inference

Jasjeet S. Sekhon

## Introduction

This website is for the distribution of "Matching" which is a R package for estimating causal effects by multivariate and propensity score matching. The package provides functions for multivariate and propensity score matching and for finding optimal balance based on a genetic search algorithm. A variety of univariate and multivariate tests to determine if balance has been obtained are also provided. These tests can also be used to determine if an experiment or quasi-experiment is balanced on baseline covariates. 

For an introduction to the package with documentation and examples,
please see "Multivariate and Propensity Score Matching Software with
Automated Balance Optimization: The Matching package for R." Journal
of Statistical Software, 42(7): 1-52. 2011.

## How to install

A version is on CRAN. The latest development version can be installed directly from Github
using [devtools](https://github.com/hadley/devtools). 

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("JasjeetSekhon/Matching")
```

The package contains compiled code, and you must have a development
environment to install the development version. (Use
`devtools::has_devel()` to check whether you do.) If no development
environment exists, Windows users download and install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and macOS
users download and install
[Xcode](https://itunes.apple.com/us/app/xcode/id497799835).
