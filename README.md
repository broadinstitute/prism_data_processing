<p align="center">
  <img src="img/BroadInstLogoforDigitalRGB.png" width="300" hspace="40"/>
  <img src="img/prism_logo_tagline_side.png" width="400" />
</p>

# Prism Data Processing

Public version of the data processing pipeline for [PRISM](https://www.theprismlab.org/) medium throughput screens (MTS). For use by collaborators to regenerate tables and plots correlating drug response to cell line features. Public cell line data for equivalent analysis is available on the [DepMap Portal](https://depmap.org/portal/).

For more general biomarker analysis see the `cdsr_biomarker` [package](https://github.com/broadinstitute/cdsr_biomarker) which is based on models in the `cdsr_models` [package](https://github.com/broadinstitute/cdsr_models).

This repository contains 3 primary scripts:

1. `make_logMFI.R`
2. `MTS_Data_Processing.R`
3. `MTS_Analysis.R`

**FIRST** run `setup.R` either in RStudio or terminal to install  packages required for analysis. For shell execute:
```bash
$ cd ./prism_data_processing
$ Rscript src/setup.R
```

For more information and FAQs see the [info](./info) folder.

### `make_logMFI.R`
Only necessary for processing raw files downloaded from [clue.io](https://clue.io/).

Converts raw .gctx and .txt files downloaded from clue.io to readable logMFI.csv. This file contains raw log2 median fluorescence intensity (MFI) data for each cell line at each treatment. Requires:
- `PR300_LMFI.gctx`: readout for PR300
- `PR500_LMFI.gctx`: readout for PR500
- `PR300_inst_info.txt`: treatment info for PR300
- `PR300_cell_info.txt`: cell line info for PR300
- `PR500_inst_info.txt`: treatment info for PR500
- `PR500_cell_info.txt`: cell line info for PR500
- `skipped_wells.csv`: file indicating which wells did not receive compound (optional)

### `MTS_Data_Processing.R`

Does pre-processing based on raw median fluorescent intensity (MFI) values and generates tables with data on QC, viabilities, and dose-response parameters, as well as figures showing dose response curves.

Steps of pre-processing outlined in [`MTS_pipeline.md`](.info/MTS_pipeline.md).


### `MTS_Analysis.R`

Generates biomarker analysis for the processed data, including, univariate and multivariate analyses. Requires a directory of expression data (RNA, mutations, etc.) and the results of `MTS_Data_Processing.R`. See below for more details on each analysis function. Relies on `analysis_functions.R`.

**NOTE:** in order to run biomarker analysis, files containing omics data for cell lines must be downloaded from [DepMap](https://depmap.org/portal/download/all/). In particular we recommend the latest versions of:
- `Achilles_gene_effect.csv` (CRISPR dependencies)
- `CCLE_expression.csv` (gene expression)
- `primary-screen-replicate-collapsed-logfold-change.csv` (repurposing)
- `D2_Achilles_gene_dep_scores.csv` (shRNA)
- `CCLE_metabolomics_20190502.csv` (metablomics)
- `CCLE_RPPA_20181003.csv` (proteomics)
- `CCLE_miRNA_20181103.gct` (miRNA)
- `CCLE_gene_cn.csv` (copy number)
- `CCLE_mutations.csv` (mutations)
- `sample_info.csv` (lineages)

Once these files are downloaded, use [`biomarker_tables.R`](src/biomarker_tables.R) to convert tables to matrix form (this may require tweaking due to changes in file structures). This R script also generates two combined datasets: `x-ccle.csv` and `x-all.csv`. These are used for multivariate models and are based on CCLE data and all DepMap data respectively.

### `MTS_functions.R` and `analysis_functions.R`

These files contain helper functions used in the scripts above and are sourced at the beginning of each (to install necessary packages and define functions). `analysis_functions.R` contains functions used in `MTS_Analysis.R` to generate biomarker analyses, while `MTS_functions.R` contains functions used in `MTS_Data_Processing.R`.

In `analysis_functions.R`, each function takes a matrix of features (X) and a vector of responses (y) as input. `multivariate` fits both an elastic net and random forest to the data.

---
_This repository is maintained by [Cancer Data Science](https://www.cancerdatascience.org/) at the Broad Institute of MIT and Harvard_
