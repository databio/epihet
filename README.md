# Epihet: Epigenetic heterogeneity from DNA methylation
[![Build Status](https://travis-ci.org/databio/RPIM.svg?branch=master)](https://travis-ci.org/databio/RPIM)

`Epihet` calculates a the PIM score, which measures the epigenetic heterogeneity in a bisulfite sequencing sample. Under the assumption that a homogeneous sample will have mostly CpGs with either 100% or 0% DNA methylation, it follows that the proportion of sites that differ from these two extremes can be used as a measure of sample heterogeneity.



### Installing

Install the development version of epihet directly from GitHub:

```{r}

devtools::install_github("databio/epihet")

```



or locally after downloading/cloning the source code:

```{r}

install.packages("path/to/epihet/directory", repos=NULL, type="source")

```

### How to use

Please see the [vignettes](/vignettes) for details.
