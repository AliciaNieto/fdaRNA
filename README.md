The *fdaRNA* package.
=====================

This package is dedicated to two central aspects of gene expression
preprocessing: normalization and outlier detection. In this package,
those methodologies are based on statistical depth. The outlier
detection methodology can be generally applied.

For more details please contact: [Alicia Nieto-Reyes](alicia.nieto@unican.es) 
and [Javier Cabrera](xavier.cabrera@gmail.com).

How to install fdaRNA?
-------------------------

The current development version can be downloaded from GitHub via

``` r
if (!requireNamespace("remotes")) install.packages("remotes")

remotes::install_bioc("affy")
remotes::install_github("asael697/fdaRNA",dependencies = TRUE)
```