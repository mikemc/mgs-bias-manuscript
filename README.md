# Data and analysis for McLaren, Willis, and Callahan (2019)

McLaren MR, Willis AD, Callahan BJ. 2019. Consistent and correctable bias in
metagenomic sequencing measurements. bioRxiv 559831.
https://www.biorxiv.org/content/10.1101/559831v1

## R package

This repository is an [R package](http://r-pkgs.had.co.nz/). It can be
installed or loaded without installation with `devtools`.  Once loaded, all
custom functions (defined in `R/`) used for our analysis will be available and
all data objects used for our main statistical analyses (stored in `data/`) can
be loaded with the `data` function.

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

The `vignettes/` folder contains a demonstration on simulated data of the bias
estimation procedure we use in our main analysis.  A rendered version can be
found [here](https://mikemc.github.io/2019-bias-manuscript/vignettes/bias-estimation-demo.html). 
