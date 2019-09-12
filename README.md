# prism\_data\_processing
Public version of the data processing pipeline for PRISM medium throughput
screens (MTS)

This repository contains 4 primary scripts:

1. `MTS_Data_Processing.R`
2. `Univariate_Continuous.R`
3. `Univariate_Discrete.R`
4. `Multivariate.R`

### [`MTS_Data_Processing.R`](./MTS_Data_Processing.R)

Does pre-processing based on raw MFI values and generates tables and dose
response curves.

Steps of pre-processing outlined in [`MTS_pipeline.md`](./MTS_pipeline.md)


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

This file contains helper functions used in the scripts above and is sourced at the beginning of each.
