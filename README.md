
MultiJGL project is about developing a general framework for the inference of dynamical multi-dependency (linear and nonlinear) networks across space and time. 
The current implementation enables to use either 1) discrete and unordered, as well as 2) continuous and ordered spatio-temproal configurations. For example, 
a typical usage in biomedical applications could be:




<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

This is an example of how to list things you need to use the software and how to install them.
To get the unreviewed developmental version of multiJGL:

```r
devtools::install_github("KontioJuho/MultiJGL")
```
### Installation

1. Get a free API Key at [https://example.com](https://example.com)
2. Clone the repo
   ```sh
   git clone https://github.com/github_username/repo_name.git
   ```
3. Install NPM packages
   ```sh
   npm install
   ```
4. Enter your API in `config.js`
   ```js
   const API_KEY = 'ENTER YOUR API';
   ```

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- USAGE EXAMPLES -->
## Acute myeloid leukemia example: 

**Joint estimation of cytogenetic risk-group specific gene networks** 

Download clinical, molecular, and  mRNA data
using the providedXenaprep function. and define the grouping variable (AMLgroups) for this analysis as an individual-level 
cytogenetic risk classification (adverse, intermerdiate, favorable). 

```r

Xenadata <- multiJGL::Xenaprep("LAML")

Xenadata$clin.molsubtype.mRNA %>% 
  rename(., cytgenrisk = acute_myeloid_leukemia_calgb_cytogenetics_risk_category) %>% 
  filter(., !is.na(cytgenrisk)) -> AMLtcga

AMLgroups <- AMLtcga$cytgenrisk

```
Identify the most frequently mutated genes in AML samples, e.g., by using _maftools_ R package
Extract the expression levels (mRNA) of these genes from the downloaded AML dataset 
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
 select(som_int$gene1) %>%
  scale(, center = TRUE, scale = TRUE) -> genes

```

```r
#It seems that five covariates are left-skewed
 multiJGL::anomaly_check(genes)

#It seems that five covariates are left-skewed
 apply(genes, 2, function(x) skewness(x)) %>% 
   filter(skewness < -1) %>%
    rownames( )  -> skew.genes
 #Transformation
  genes[,skew.genes] <- log10(max(genes[,skew.genes]+1) - genes[,skew.genes])
  
  #Recheck
  anomaly_check(genes)

  
```

```r
networks <- multiJGL(node.covariates = genes,
                           grouping.factor = AMLgroups,
                           lin_lambda1 = 0.1, lin_lambda2 = 0.025,
                           nonlin_lambda1 = 0.1,nonlin_lambda2 = 0.025)
```


_For more examples, please refer to the [Documentation](https://example.com)_

<p align="right">(<a href="#top">back to top</a>)</p>
