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
generate_geno(ancestries=c("EUR","EAS"), #ancestries
              replication=3, #number of replications,
              chr_set=20, #simulated chrmosomes
              sample_size=10000,
              #qulaity control parameters
              qc_maf=0.05,
              qc_hwe=1e-7,
              qc_geno=0.05,
              qc_mind=0.05,
              ouput_path=getwd()
              )
```

#### Feature 2. Generate omics data

``` R
library(simGO)
geno2omics(h2_omics=rep(0.1,367),  #h2 of omics contributed by geno 
           cov_omics=InterSIM::cov.M, # covariance across omics data
           mean.omics.cluster=InterSIM::mean.M,
           n.cluster=4, #cluster of the number of data
           p.cluster=c(0.25,0.25,0.25,0.25), #proportions of the sample within each cluster
           delta.omics.cluster=5,
           align_block_index,
           p.deomics=0.2, #proportions of the omics features are differentialy expressed
           geno=genotype_matrix, #omics1 data, should have been scaled 
           seed=1998
                  )                 
```

#### Feature 3. Generate phenotype data

``` R
library(simGO)
#Sparse
prefix_bin="/restricted/projectnb/adiposity/qsmei/SimGO/Genotype/EUR/EUR" #example of genotype data
data1=simGO(prefix_bin=prefix_bin,nTraits=2,
            h2SNPs=list(0.3,0.5),  #heratability accounted by genotype data 
            SNPs_archi="Sparse", #Infinite,Sparse,BSLMM,BayesR
            SNPs_causal_var=1,#Variances of SNPs within each group,
            pSNPs_causal_var=0.01,#proportions of SNPs with specified variance within each group,correspond to SNPs_var, 
            #nSNPs_causal_var=1000,#number of SNPs with specified variance within each group,correspond to SNPs_var, 
            #if input of SNPs(pSNPs_causal_var <1 or nSNPs_causal_var < nSNPs),the rest SNPs are stand for zero-effect SNPs 
            SNPs_causal_mode="Random", #Random, means random sampling, ignore LD
            sSNPs_negative=-1, #{√(2pq)}^s, control the degree of negative selection,
            #pSNPs_causal_shared=0.01,
            #nSNPs_causal_shared=1000, # number of shared causal SNPs across all causal SNPs
            SNPs_cor=0.8,
            Res_cor=0.1,
            seed=1998,
            output_prefix = "Sparse_0.3",
            output_path=getwd() # default is the current working path
            )

```
