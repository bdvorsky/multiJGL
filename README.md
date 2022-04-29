


<!-- GETTING STARTED -->
## Getting Started


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
