# simGO <img src="https://img.shields.io/badge/Issues-%2B-brightgreen.svg" /><img src="https://img.shields.io/badge/license-GPL3.0-blue.svg" />
### A Comprehensive Simulation Framework for Modeling Relationships Among Genotypes, Omics, and Phenotypes 
 <img src="https://github.com/qsmei/simGO/blob/main/Supplementary/Workflow.jpg" alt="logo-simGO" align="center"/>  
## Contents

-   [OVERVIEW](#overview)

-   [GETTING STARTED](#getting-started)

    -   [Installation](#installation)
    -   [Features](#features)

-   [USAGE](#usage)

------------------------------------------------------------------------

### OVERVIEW 
`simGO` is comprehensive Simulation Framework for Modeling Relationships Among Genotypes, Omics, and Phenotypes

If you have suggestion or question, please contact: quanshun1994@gmail.com ! 

## GETTING STARTED

### 🙊Installation

`simGO` links to R packages `Rcpp`, `RcppArmadillo` , `data.table` and  `bigmemory` .  These dependencies should be installed before installing `simGO`.  

```R
install.packages(c("Rcpp", "RcppArmadillo","RcppProgress","data.table","bigmemory","R6"))
if (!require(devtools)) install.packages("devtools")#if devtools not already installed
```
**👉 Note: In the analysis of DMU  and BLUPF90 , we need to download software DMU  ([DMU download website](https://dmu.ghpc.au.dk/dmu/))  and BLUPF90 previously ([BLUPF90 download website](http://nce.ads.uga.edu/html/projects/programs/)). For convenience, we have encapsulated  the basic module of DMU and BLUPF90 in package `blupADC`.**  

 **For commercial use of DMU and BLUPF90,  user must contact the author of DMU and BLUPF90 !!!** 

#### Install the latest version of simGO
```R
remotes::install_github("qsmei/simGO")
```
or
```R
devtools::install_github("qsmei/simGO")
```

After installed successfully, the `simGO` package can be loaded by typing

``` {.r}
library(simGO)
```

**Note**: In terms of the relationship matrix construction, we highly recommend Microsoft R Open(faster than traditional R many times)

### 🙊Features

-   Feature 1. Generate genotype data
-   Feature 2. Generate omics data
-   Feature 3. Generate phenotype data


## Usage

``` {.r}
system.file("extdata", package = "simGO") # path of provided files
```

#### Feature 1. Generate genotype data
``` R
library(simGO)


```

#### Feature 2. Generate omics data

``` R
library(simGO)                  
```

#### Feature 3. Generate phenotype data

``` R
library(simGO)
```
