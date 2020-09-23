# Ovarian cancer subtype heterogeneity

This repository contains code and analysis scripts that were used for
conducting analysis and producing the figures of the publication:

Ludwig Geistlinger, Sehyun Oh, Marcel Ramos, Lucas Schiffer, Rebecca S. LaRue, 
Christine M. Henzler, Sarah A. Munro, Claire Daughters, Andrew C. Nelson, 
Boris J. Winterhoff, Zenas Chang, Shobhana Talukdar, Mihir Shetty, Sally A.
Mullany, Martin Morgan, Giovanni Parmigiani, Michael J. Birrer, Li-Xuan Qin,
Markus Riester, Timothy K. Starr and Levi Waldron.
[Multi-omic analysis of subtype evolution and heterogeneity in high-grade serous ovarian carcinoma](https://doi.org10.1158/0008-5472.CAN-20-0521).
Cancer Research, 2020. doi: 10.1158/0008-5472.CAN-20-0521.

Analysis as documented in this repository was carried out using
[R-4.0.2](https://cran.r-project.org/) and
[Bioconductor-3.11](https://www.bioconductor.org/install/).

**Installation:**

Assuming R and Bioconductor are installed:

1. Download the code from GitHub using the `Download ZIP` option.
2. Unzip the downloaded `subtypeHeterogeneity-master.zip`. 
3. Rename the unzipped directory to `subtypeHeterogeneity`.
4. Build the package eg. via `R CMD build subtypeHeterogeneity`.
5. Install the package via `R CMD INSTALL subtypeHeterogeneity_1.0.0.tar.gz`. 

## Absolute somatic copy number alteration analysis (Figures 2 and 3)

The code for reproducing Figures 2 and 3 can be found in the
[vignettes/subtypeHeterogeneity.Rmd](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/vignettes/subtypeHeterogeneity.Rmd) file.

## 10x Genomics single-cell RNA-seq analysis (Figures 4 and 5) 

The code for reproducing Figures 4 and 5 can be found in the
[vignettes/singlecell.Rmd](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/vignettes/singlecell.Rmd) file.


