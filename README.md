<!--
---
output:
  html_document:
    keep_md: yes
---

-->

# *treeWAS*: A phylogenetic tree-based approach to genome-wide association studies in microbes


<!-- ########################################################################################################## -->
## Introduction
<!-- ########################################################################################################## -->

The *treeWAS* R package allows users to apply our phylogenetic tree-based appraoch to Genome-Wide Association Studies (GWAS) to microbial genetic and phenotypic data. In short, *treeWAS* aims to identify genetic variables that are statistically associated with a phenotypic trait, while correcting for the confounding effects of population structure and recombination. *treeWAS* is applicable to both bacterial and viral genetic data and to both binary and continuous phenotypes. The approach adopted within *treeWAS* is described fully in our paper, available on [bioRXiv](http://www.biorxiv.org/content/early/2017/05/22/140798).



***

<!-- ########################################################################################################## -->
## Installation
<!-- ########################################################################################################## -->

*treeWAS* is currently hosted on GitHub at <https://github.com/caitiecollins/treeWAS>.  
<!-- ([https://github.com/caitiecollins/treeWAS](https://github.com/caitiecollins/treeWAS)).-->

The most up-to-date version of *treeWAS* can be easily installed directly within R, using the `devtools` package: 


```{r, eval=FALSE, highlight=TRUE}
## install devtools, if necessary:
install.packages("devtools", dep=TRUE)
library(devtools)

## install treeWAS from github:
install_github("caitiecollins/treeWAS/pkg", build_vignettes = TRUE)
library(treeWAS)
```


***

<!-- ########################################################################################################## -->
## Documentation
<!-- ########################################################################################################## -->

Documentation on how to use *treeWAS* can be found on GitHub in [the Wiki](https://github.com/caitiecollins/treeWAS/wiki). 


The Wiki contains sections on [The Method](https://github.com/caitiecollins/treeWAS/wiki/1.-How-treeWAS-Works) behind *treeWAS*, 
the [Data & Data Cleaning](https://github.com/caitiecollins/treeWAS/wiki/2.-Data-&-Data-Cleaning) required, 
the [treeWAS Function & Arguments](https://github.com/caitiecollins/treeWAS/wiki/3.-treeWAS-Function-&-Arguments), 
a guide to [Interpreting Output](https://github.com/caitiecollins/treeWAS/wiki/4.-Interpreting-Output) returned by *treeWAS*, 
functions to facilitate [Integration with ClonalFrameML](https://github.com/caitiecollins/treeWAS/wiki/5.-ClonalFrameML-Integration), 
and information describing how to flag [Bugs & Features](https://github.com/caitiecollins/treeWAS/wiki/6.-Bugs-&-Features).



Once you have installed and loaded the *treeWAS* package, you can also find this information in the vignette. 
To open the vignette from within R (recommended if any formatted elements are not rendering properly in the wiki), 
run `browseVignettes` and click on the `HTML` hyperlink:


```{r, eval=FALSE}
browseVignettes("treeWAS")
```

***
