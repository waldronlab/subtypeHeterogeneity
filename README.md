# Ovarian cancer subtype heterogeneity

This repository contains code and analysis scripts that were used for
conducting analysis and producing the figures of the publication:

Ludwig Geistlinger, Sehyun Oh, Marcel Ramos, Lucas Schiffer, Rebecca S. LaRue, 
Christine M. Henzler, Sarah A. Munro, Claire Daughters, Andrew C. Nelson, 
Boris J. Winterhoff, Zenas Chang, Shobhana Talukdar, Mihir Shetty, Sally A.
Mullany, Martin Morgan, Giovanni Parmigiani, Michael J. Birrer, Li-Xuan Qin,
Markus Riester, Timothy K. Starr and Levi Waldron.
[Multi-omic analysis of subtype evolution and heterogeneity in high-grade serous ovarian carcinoma](https://doi.org/10.1158/0008-5472.CAN-20-0521).
Cancer Research, 2020. doi: 10.1158/0008-5472.CAN-20-0521.

Analysis as documented in this repository was carried out using
[R-4.0.2](https://cran.r-project.org/) and
[Bioconductor-3.11](https://www.bioconductor.org/install/).

**Installation:**

Option 1: Docker

[Docker image](https://hub.docker.com/repository/docker/ludwigg/enrichomics) 
that comes with all dependencies already installed

Option 2: Local installation

Assuming [R](https://cran.r-project.org/) and 
[Bioconductor](https://www.bioconductor.org/install/) 
are installed:

```
BiocManager::install("waldronlab/subtypeHeterogeneity",
                     dependencies = TRUE,
                     build_vignettes = TRUE)
```

## Absolute somatic copy number alteration analysis (Figures 2 and 3)

The code for reproducing Figures 2 and 3 can be found in the
[vignettes/subtypeHeterogeneity.Rmd](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/vignettes/subtypeHeterogeneity.Rmd) file.

To view the fully rendered HTML document containing literate programming output
see either [here](https://waldronlab.io/subtypeHeterogeneity/articles/subtypeHeterogeneity.html), 
or from within your local R installation:

```
browseVignettes("subtypeHeterogeneity")
```

and select the `Absolute SCNA analysis` HTML vignette.

## 10x Genomics single-cell RNA-seq analysis (Figures 4 and 5) 

The code for reproducing Figures 4 and 5 can be found in the
[vignettes/singlecell.Rmd](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/vignettes/singlecell.Rmd) file.

To view the fully rendered HTML document containing literate programming output
from within R:

```
browseVignettes("subtypeHeterogeneity")
```

and select the `Single cell analysis` HTML vignette.

