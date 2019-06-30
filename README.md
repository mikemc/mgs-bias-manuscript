# Data and analysis for McLaren, Willis, and Callahan (2019)

McLaren MR, Willis AD, Callahan BJ. 2019. Consistent and correctable bias in
metagenomic sequencing measurements. bioRxiv 559831.
https://www.biorxiv.org/content/10.1101/559831v1

## This repository

This repository contains the code and data for reproducing the analysis in our
manuscript. It is structured as an [R package](http://r-pkgs.had.co.nz/) in
order to better manage and document the custom R functions we use in our
analysis; see a [discussion of this
approach](https://github.com/ropensci/rrrpkg). The primary goal of this
repository is thus to enable others to understand and reproduce our results,
and not to provide an R package to be used directly in other applications.
Once the package is loaded, all custom functions (defined in `R/`) used for our
analysis will be available and all data objects used for our main statistical
analyses (stored in `data/`) can be loaded with the `data` function. 

To reproduce our analysis, download the package from Github or with
```
git clone https://github.com/mikemc/2019-bias-manuscript
```
and `knit` or run the R-markdown documents in `analysis/` (see below).  The
package can be loaded without installing it to your R library using the
function `devtools::load_all`. This is what is done by the Rmd documents in
`analysis/`. Alternatively, you can install it and load it like a normal R
package:
```r
# Load package after downloading (recommended)
devtools::load_all("/path/to/local/2019-bias-manuscript")

# Alternatively, install the package with
# devtools::install_git("/path/to/local/2019-bias-manuscript")
# or
# devtools::install_github("mikemc/2019-bias-manuscript")
# and load as a normal R package
# library(BiasManuscript)

data("brooks2015_sample_data")
data("brooks2015_counts")
data("brooks2015_species_info")

library(dplyr)

brooks2015_sample_data %>%
    group_by(Mixture_type) %>%
    summarize(Num_samples = n())
#> # A tibble: 3 x 2
#>   Mixture_type Num_samples
#>   <chr>              <int>
#> 1 Cells                 80
#> 2 DNA                   80
#> 3 PCR_product           80
```

## Data

The scripts we used to download and/or generate the necessary sample metadata,
16S and metagenomic taxonomic profiles, and taxon information for our analyses
are in `data-raw/`.  This folder also contains scripts that clean the data and
save it as `.rda` (R data) objects that can be loaded with the `data` function
once the R package is loaded.  These objects serve as the starting point for
subsequent analyses.

## Analysis

Analyses are contained in R-markdown documents in `analysis/`. Versions
rendered to html can be seen at
* [Analysis of the Brooks et al (2015) dataset](https://mikemc.github.io/2019-bias-manuscript/analysis/brooks2015-analysis.html)
* [Analysis of the Costea et al (2017) dataset](https://mikemc.github.io/2019-bias-manuscript/analysis/costea2017-analysis.html)
* [Calculations for conceptual examples discussed in the manuscript](https://mikemc.github.io/2019-bias-manuscript/analysis/conceptual-examples.html)
* [Estimating genome statistics for the Brooks et al (2015) species](https://mikemc.github.io/2019-bias-manuscript/analysis/brooks2015-species-info.html)
