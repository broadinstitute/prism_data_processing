<p align="center">
  <img src="BroadInstLogoforDigitalRGB.png" width="300" hspace="40"/>
  <img src="prism_logo_tagline_side.png" width="400" />
</p>

# Prism Data Processing

Public version of the data processing pipeline for PRISM medium throughput screens (MTS). For use by collaborators to regenerate tables and plots correlating drug response to cell line features. Public cell line data for equivalent analysis is available on the [DepMap Portal](https://depmap.org/portal/).

This repository contains 5 primary scripts:

1. `MTS_Data_Processing.R`
2. `MTS_Analysis.R`
3. `Univariate_Continuous.R`
4. `Univariate_Discrete.R`
5. `Multivariate.R`

**FIRST** run [`requirements.R`](./requirements.R) either in RStudio or terminal to install  packages required for analysis. For shell execute:
```bash
$ Rscript requirements.R
```

### [`MTS_Data_Processing.R`](./MTS_Data_Processing.R)

Does pre-processing based on raw median fluorescent intensity (MFI) values and generates tables with data on QC, viabilities, and dose-response parameters, as well as figures showing dose response curves.

Steps of pre-processing outlined in [`MTS_pipeline.md`](./MTS_pipeline.md)


### [`MTS_Analysis.R`](./MTS_Analysis.R)

Generates biomarker analysis for the processed data, including, univariate and multivariate analyses. Requires a directory of expression data (RNA, mutations, etc.) and the results of `MTS_Data_Processing.R`. See below for more details on each analysis function. Relies on `analysis_functions.R`.

Note: it is recommended to use this script over individual analyses as it has the most up to date methods. To run a single analysis (e.g. correlation of a feature like dependency score with an assay result like AUC), use the `correlate`, `discrete_test`, or `multivariate` functions in `analysis_functions.R`.

### [`MTS_functions.R`](./MTS_functions.R) and [`analysis_functions.R`](./analysis_functions.R)

These files contain helper functions used in the scripts above and are sourced at the beginning of each (to install necessary packages and define functions). `analysis_functions.R` contains functions used in `MTS_Analysis.R` to generate biomarker analyses, while `MTS_functions.R` contains functions used in `MTS_Data_Processing.R`.

In `analysis_functions.R`, each function takes a matrix of features (X) and a vector of responses (y) as input. `multivariate` fits both an elastic net and random forest to the data.

---
_This repository is maintained by [Cancer Data Science](https://www.cancerdatascience.org/) at the Broad Institute of MIT and Harvard_
