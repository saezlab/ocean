# ocean
R package for metabolic enzyme enrichment anaylsis

![alt text](https://github.com/saezlab/ocean/blob/master/ocean_logo.001.png?raw=true)

## Tutorial

Instal the package (from github with devtools) :

```r
## If needed instal devtool package
install.packages("devtools")

## instal COSMOS
library(devtools)
install_github("saezlab/ocean")
```

You can then run the tutorial scripts with a kidney cancer toy metabolomic dataset: https://github.com/saezlab/ocean/blob/master/tutorial_submodules.R

Ocean metabolic model is based on a curated reduced model of human metabolism, see: Masid M, Ataman M & Hatzimanikatis V (2020) Analysis of human metabolism by reducing the complexity of the genome-scale models using redHUMAN. Nat Commun 11: 2821 https://www.nature.com/articles/s41467-020-16549-2#Tab1
