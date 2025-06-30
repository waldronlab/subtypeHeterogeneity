# Ovarian cancer subtype heterogeneity

## How can I obtain cell metadata inclluding cell type and subtype annotation (Figure 4)?

**Option 1: CSV files**

Available as CSV files in the `inst/extdata` directory of this repository: 

* [T59](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/inst/extdata/T59_cell_metadata.csv)
[T76](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/inst/extdata/T76_cell_metadata.csv)
[T77](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/inst/extdata/T77_cell_metadata.csv)
[T89](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/inst/extdata/T89_cell_metadata.csv)
[T90](https://github.com/waldronlab/subtypeHeterogeneity/blob/master/inst/extdata/T90_cell_metadata.csv)

**Option 2: R/Bioconductor SingleCellExperiment objects**

Assuming familiarity with R and with the [SingleCellExperiment](https://bioconductor.org/books/release/OSCA.intro/the-singlecellexperiment-class.html) container, all you would need to do is to carry out lines 290-310 here:

https://github.com/waldronlab/subtypeHeterogeneity/blob/0f275d2f5a582046a1f36f4150970157da6e6fe1/vignettes/singlecell.Rmd#L290

and access the colData of the obtained SingleCellExperiments which includes, among other things, the cell type annotation.

For example here for Tumor T59:

```R
> colData(sces[["T59"]])[,c("Barcode", "subtype", "celltype")]
DataFrame with 13040 rows and 3 columns
                 Barcode     subtype    celltype
             <character> <character> <character>
1     AAACCTGAGCTGCCCA-1         DIF         MYE
2     AAACCTGAGTCATCCA-1         DIF         MYE
3     AAACCTGCAAGCCCAC-1         IMR         MYE
4     AAACCTGCAAGCGCTC-1         DIF         EPI
5     AAACCTGCACGTAAGG-1         DIF         EPI
...                  ...         ...         ...
13036 TTTGTCATCCCTTGCA-1         DIF         MYE
13037 TTTGTCATCCGTACAA-1         DIF         MYE
13038 TTTGTCATCCTCCTAG-1         DIF         EPI
13039 TTTGTCATCGGTCTAA-1         DIF         EPI
13040 TTTGTCATCTAACTCT-1         DIF         MYE
```


## Overview

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

## Installation:

**Option 1: Docker**

[Docker image](https://hub.docker.com/repository/docker/ludwigg/subtypeheterogeneity) 
that comes with all dependencies already installed

**Option 2: Local installation**

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
see either [here](https://waldronlab.io/subtypeHeterogeneity/articles/singlecell.html), 
or from within your local R installation:

```
browseVignettes("subtypeHeterogeneity")
```

and select the `Single cell analysis` HTML vignette.

