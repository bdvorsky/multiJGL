# MultiJGL <img src="man/logo/multijgl.png" align="right" height = 150/>

<!-- badges: start -->
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

## Overview

MultiJGL project is about developing a general framework for the inference of multi-dependency (linear and nonlinear) networks across multiple sample groups. The current implementation introduces a nonlinear extension of the popular joint graphical lasso (JGL) enabling the identifcation of both linear and nonlinear dependencies between nodes in a hierarchical manner. 

#### Output/Result: 

Description of input data and description of algorithm + output according to the picture bellow:

Here comes your overview picture:

<img src="example/new-figure-10-1024x576.png"  height = 300/>


<!-- GETTING STARTED -->
## Install the package
To install the developer version (status: unreviewed) from github:
```r
library(devtools)
devtools::install_github("KontioJuho/multiJGL")
``` 

## Cheat sheet
Is it viable here?

<a href="https://github.com/rstudio/cheatsheets/blob/main/data-import.pdf"><img src="https://raw.githubusercontent.com/rstudio/cheatsheets/main/pngs/thumbnails/data-import-cheatsheet-thumbs.png" width="630" height="252"/></a>

<!-- USAGE EXAMPLES -->
## Usage: 

MultiJGL creates lorem ipsum lorem ipsum lorem ipsum lorem ipsum lorem ipsum

Please follow these documented use cases for multiJGL for undestanding how to use this package:

- [Use case 1]()
- [Use case 2]()
- [Use case 3]()

Or test the basic functionality directly with following example.


## Acute myeloid leukemia example: 

**Joint estimation of cytogenetic risk-group specific gene networks** 

Download clinical, molecular, and  mRNA data
using the providedXenaprep function. and define the grouping variable (AMLgroups) for this analysis as an individual-level 
cytogenetic risk classification (adverse, intermerdiate, favorable). 

```r

library('multiJGL')
Xenadata <- multiJGL::Xenaprep("LAML")

Xenadata$clin.molsubtype.mRNA %>% 
  rename(., cytgenrisk = acute_myeloid_leukemia_calgb_cytogenetics_risk_category) %>% 
  filter(., !is.na(cytgenrisk)) -> AMLtcga

AMLgroups <- AMLtcga$cytgenrisk

```
Identify the most frequently mutated genes in AML samples, e.g., by using _maftools_ R package, and
extract the expression levels (mRNA) of these genes from the downloaded AML dataset: 
```r

library(maftools)
library(e1071)
#Set path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')

#Query the names of 15 most frequently mutated genes in AML
read.maf(maf = laml.maf, clinicalData = laml.clin) %>% 
  somaticInteractions(maf = , top = 15) -> som_int 

AMLtcga %>%
  dplyr::select(.,som_int$gene1) %>%
  scale(., center = TRUE, scale = TRUE) -> genes

```
Validate the normality assumption with a built-in _data_check_ function and run an appropriate transformation function if needed:
```r
#It seems that five covariates are left-skewed
multiJGL::data_check(genes)

library(e1071)
#It seems that five covariates are left-skewed
apply(genes, 2, function(x) skewness(x)) %>% 
  as.data.frame() %>% 
  filter(. < -1) %>%
  rownames( )  -> skew.genes
#Transformation
genes[,skew.genes] <- scale(log10(max(genes[,skew.genes]+1) - genes[,skew.genes]))

#Recheck
multiJGL::data_check(genes)

  
```
Estimate both linear and nonlinear network structures with _multiJGL_ function (with the default parameter values):
```r
networks <- multiJGL(node.covariates = genes,
                           grouping.factor = AMLgroups,
                           lin_lambda1 = 0.1, lin_lambda2 = 0.025,
                           nonlin_lambda1 = 0.1,nonlin_lambda2 = 0.025)
```


_For more examples, please refer to the [Documentation]()_


## Getting help: 

- contact the developer
- ask question...

<p align="right">(<a href="#top">back to top</a>)</p>
