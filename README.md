# prism\_data\_processing
Public version of the data processing pipeline for PRISM medium throughput
screens (MTS)

This repository contains 5 primary scripts:

1. `MTS_Data_Processing.R`
2. `MTS_Analysis.R`
3. `Univariate_Continuous.R`
4. `Univariate_Discrete.R`
5. `Multivariate.R`

**FIRST** open R Studio and run [`requirements.R`](./requirements.R) to install
required packages.

### [`MTS_Data_Processing.R`](./MTS_Data_Processing.R)

Does pre-processing based on raw MFI values and generates tables and dose
response curves.

Steps of pre-processing outlined in [`MTS_pipeline.md`](./MTS_pipeline.md)


### [`MTS_Analysis.R`](./MTS_Analysis.R)

Generates biomarker analysis for the processed data, including, univariate and
multivariate analyses. Requires a directory of expression data (RNA, mutations,
etc.) and the results of `MTS_Data_Processing.R`. See below for more details on
each analysis function. Relies on `analysis_functions.R`.

Note: it is recommended to use this script over individual analyses as it has
the most up to date methods.

### [`Univariate_Continuous.R`](./Univariate_Continuous.R)

Calculates univariate correlations between a given continuous feature of a cell
line (e.g. gene expression) and response to a given compound (any response
metric such as IC50, AUC, or LMFI can be used).


### [`Univariate_Discrete.R`](./Univariate_Discrete.R)

Runs t-tests for groups of cell lines defined by discrete variables (e.g.
lineage), comparing responses between groups.


### [`Multivariate.R`](./Multivariate.R)

Trains a random forest and elastic net models on the data given a feature and a
response. These models then output the relative importance of each feature in
predicting response.

### [`MTS_functions.R`](./MTS_functions.R)

This file contains helper functions used in the scripts above and is sourced at
the beginning of each.
