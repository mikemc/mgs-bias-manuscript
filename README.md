# Data and analysis for McLaren, Willis, and Callahan (2019)

McLaren MR, Willis AD, Callahan BJ. 2019. Consistent and correctable bias in
metagenomic sequencing measurements. bioRxiv 559831.
https://www.biorxiv.org/content/10.1101/559831v2

## This repository

This repository contains the code and data for reproducing the analysis in our
manuscript. It is structured as an [R package](http://r-pkgs.had.co.nz/), as 
explained [here](https://github.com/ropensci/rrrpkg). To reproduce our 
analysis, first install the manuscript version of the 
[metacal R package](https://github.com/mikemc/metacal)
```r
# install.packages("devtools")
devtools::install_github("mikemc/metacal@v0.1.0-manuscript")
```
Then, download this package from GitHub or by running
```
git clone https://github.com/mikemc/mgs-bias-manuscript
```
You can then `knit` or run the R-markdown documents in `analysis/`, which are
described below. These documents include code to load this package with
`devtools::load_all()`, so you do not need to install this package itself.
Various other R packages are needed to run the code in the `analysis/`
documents; these are listed in the "Imports" field of the `DESCRIPTION` file
and can be installed all at once with
```r
devtools::install_deps("path/to/mgs-bias-manucript")
```

## Data

The scripts we used to download and/or generate the necessary sample metadata,
16S and metagenomic taxonomic profiles, and taxon information for our analyses
are in `data-raw/`. This folder also contains scripts that clean the
data and save it as `.rda` (R data) objects that can be loaded with the 
`data()` function once the R package is loaded; these objects serve as the 
starting point for subsequent analyses. An explanation of how to use these 
scripts is given in the directory's Readme file.

## Analysis

Analyses are contained in R-markdown documents in `analysis/`. 
Versions already rendered to html can be seen at
* [Calculations for conceptual examples discussed in the manuscript](https://mikemc.github.io/mgs-bias-manuscript/analysis/conceptual-examples.html)
* [Estimating genome statistics for the Brooks et al (2015) species](https://mikemc.github.io/mgs-bias-manuscript/analysis/brooks2015-species-info.html)
* [Analysis of the Brooks et al (2015) dataset](https://mikemc.github.io/mgs-bias-manuscript/analysis/brooks2015-analysis.html)
* [Analysis of the Costea et al (2017) dataset](https://mikemc.github.io/mgs-bias-manuscript/analysis/costea2017-analysis.html)

