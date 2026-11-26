# conformeR

Conformal clustering of gene differential expression.

`conformeR` is a wrapper around the [LEMUR](https://github.com/const-ae/lemur) 
R package that adds
uncertainty quantification to LEMUR's neighborhood on predicted differential 
expression, using conformal prediction, and without
altering the underlying model.

<img src="man/figures/logo.png" align="right" height="139" alt="" />








## Installation

`conformeR` is available from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("conformeR")
```

The development version can be installed from GitHub:

```r
BiocManager::install("juslecl/conformeR")
```

## Usage

```r
library(conformeR)
```

See `?conformeR` for the main entry point.

## Vignette
Open it with:
 
```r
browseVignettes("conformeR")
```
 
or, from the R console:
 
```r
vignette("conformeR")
```

## Paper
Stay tuned!
