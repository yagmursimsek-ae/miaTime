# Microbiome time series manipulation with miaTime <img src="man/figures/mia_logo.png" align="right" width="120" />

<!-- badges: start -->

[![Platforms](http://bioconductor.org/shields/availability/release/miaTime.svg)](https://bioconductor.org/packages/release/bioc/html/miaTime.html)
[![rworkflows](https://github.com/microbiome/miaTime/actions/workflows/rworkflows.yml/badge.svg?branch=devel)](https://github.com/microbiome/miaTime/actions)
[![Bioc-release](http://bioconductor.org/shields/build/release/bioc/miaTime.svg)](http://bioconductor.org/packages/release/bioc/html/miaTime.html)
[![Bioc-age](http://bioconductor.org/shields/years-in-bioc/miaTime.svg)](https://bioconductor.org/packages/release/bioc/html/miaTime.html)
[![Codecov test coverage](https://codecov.io/gh/microbiome/miaTime/branch/devel/graph/badge.svg)](https://codecov.io/gh/microbiome/miaTime?branch=devel)
[![Dependencies](http://bioconductor.org//shields/dependencies/release/miaTime.svg)](https://bioconductor.org/packages/release/bioc/html/miaTime.html)

<!-- badges: end -->

## Using the package

This R package can be used to analyse time series data for microbial
communities. The package is part of [miaverse](https://microbiome.github.io/), 
and is based on the `TreeSummarizedExperiment` data container.

See the [package homepage](https://microbiome.github.io/miaTime) for
example workflows.

## Installation

### Bioc-release

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("miaTime")
```

### Bioc-devel

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("miaTime")
```

### GitHub

```
remotes::install_github("microbiome/miaTime")
```

## Contributions and acknowledgments

You can find us online from [Gitter](https://gitter.im/microbiome/miaverse).

Contributions are very welcome through issues and pull requests at the
[development site](https://github.com/microbiome/miaTime). We follow a git
flow kind of approach. Development version should be done against the
`main` branch and then merged to `release` for release.
(https://guides.github.com/introduction/flow/)

**Kindly cite this work**. For citation details, see R command

```r
citation("miaTime")  
```

## Code of conduct

The project is released with a
[Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
