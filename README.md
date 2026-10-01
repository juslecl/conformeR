# conformeR

Conformalized differential expression analysis of multi-condition single-cell data.

`conformeR` is a wrapper around the [LEMUR](https://github.com/const-ae/lemur) R package that adds
uncertainty quantification to LEMUR's neighborhood on predicted differential expression, using conformal prediction, and without
altering the underlying model.

## Installation

`conformeR` is not yet on CRAN or Bioconductor. Install the development version directly from
GitHub.

```r
install.packages("remotes")
remotes::install_github("juslecl/conformeR", build_vignettes = TRUE)
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
