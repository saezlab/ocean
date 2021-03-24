# ocean
R package for metabolic enzyme enrichment anaylsis

![alt text](https://github.com/saezlab/ocean/blob/master/ocean_logo.001.png?raw=true)

## Overview

The functional insights that metabolomic data sets contain currently lies under-exploited. This is in part due to the complexity of metabolic reaction networks and the indirect relationship between reaction fluxes and metabolite abundance. Yet, footprint-based methods have been available for decades in the context of other omic data sets such as transcriptomic and phosphoproteomic. Here, we present ocEAn, a method that defines metabolic enzyme footprint from a curated reduced version of the recon2 reaction network and use them to explore coordinated deregulations of metabolite abundances with respect to their position relative to metabolic enzymes in the same manner as Kinase-substrate and TF-targets Enrichment analysis. This picture diplays the TCA cycle with deregulated metabolites and estimated metabolic enzyme deregulations in kidney cancer. You can generate a dynamic visualisation of this network and other pathways by following our tutorial: https://github.com/saezlab/ocean/blob/master/tutorial_submodules.R

![alt text](https://github.com/saezlab/ocean/blob/master/TCA_shot.png?raw=true)

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

## Citations

The manuscript associated with ocEAn will be available soon. Stay tuned !

The underlying metabolic model is based on a curated reduced model of human metabolism, see: Masid M, Ataman M & Hatzimanikatis V (2020) Analysis of human metabolism by reducing the complexity of the genome-scale models using redHUMAN. Nat Commun 11: 2821 https://www.nature.com/articles/s41467-020-16549-2#Tab1
