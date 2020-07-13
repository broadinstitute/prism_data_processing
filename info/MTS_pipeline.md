# MTS Processing Pipeline

This document outlines the steps executed by `MTS_Data_Processing.R`, an R script to process Level 2 MTS data. That script can be used to regenerate any data obtained from PRISM based on logMFI values. This information is also available in the "Data Processing" section of the report.

- **Input** is mean-fluorescence intensity of cell lines treated with small molecules and controls.
- **Output** is files describing dose-response behavior of each cell line to each drug, median-fold change of cell lines in response to small molecules, and statistical significance metrics.
- **To execute:** change `data_path` variable in `MTS_Data_Processing.R` to match local path of unprocessed data (if data is in .gctx form, use `make_logMFI.R` first).

The "resulting data" section of each step is meant to guide the user through the code and does not necessarily reflect output or files generated (only internal R variables).

### 1) Load Data
1. Reads in a .csv of logMFI (median fluorescence intensity, a proxy for cell count) data for each sample (each sample being 1 cell line treated with 1 drug)
2. Splits data into PR300 and PR500 datasets (PR300 and PR500 refer to unique sets of PRISM cell lines, referred to as cultures, which have slightly different properties)
3. Creates reference tables containing only control barcodes (pool ID is "CTLBC", these could be either negative, DMSO treated, or positive, treated with 10µM Bortezomib)
4. **Resulting data:** `PR300`, `PR500`, `PR300_barcodes`, `PR500_barcodes`

### 2) Normalize control vehicle samples (negative controls)
1. Calculates the median logMFI of each well across replicates for negative controls (DMSO treated)
2. Takes the median of calculated median to generate overall median for negative controls
3. Subtracts the well median from the sample value and adds back overall DMSO median (result is normalized samples with respect to overall DMSO median)
4. **Resulting data:** `PR300_normalized` (just negative controls), `PR500_normalized` (just negative controls)

### 3) Normalize all samples using controls
1. Fits a curve to the normalized control data relating logMFI to normalized logMFI, this defines a "normalization function" that can be applied to treatment data
2. Uses the equation to normalize non-control data
3. Recombines the newly normalized data with original data, replacing logMFI with normalized values (from here on logMFI refers to this normalized value)
4. **Resulting data:** `PR300_normalized` (all data), `PR500_normalized` (all data)

### 4) Generate statistics table for controls
1. For each control treatment pool:
 1. Calculates the median logMFI (MD)
 2. Calculates the median absolute deviation: $MAD = med(|logMFI - med(logMFI)|)$
 3. Calculates the strictly standardized mean difference: $SSMD = \frac{med_{neg} - med_{pos}}{MAD_{neg}^2 + MAD_{pos}^2}$
 4. Calculates the null-normalized mean difference: $NNMD = \frac{med_{neg} - med_{pos}}{MAD_{neg}}
 5. Calculate the error rate of logMFI distributions: error rate is determined by the accuracy of a perfect threshold classifier of control data
2. Remove samples with an error rate > 0.05 from further analysis
3. **Resulting data:** `SSMD_TABLE`

### 5) Calculate log-fold changes (LFC) with respect to negative controls
1. Only use samples with error rate < 0.05
2. Wells must be present or have a low enough error rate on at least 2 replicates
3. $LFC_{sample} = logMFI_{sample} - median(logMFI_{controls})$
4. **Resulting data:** `LFC_TABLE`

### 6) Correct for pool effects using ComBat
1. Generates a corrected LFC for each sample (`LFC.cb`) and stores in new table (separate from uncorrected)
2. **Resulting data:** `LFC_TABLE` (with new LFC.cb value)

### 7) Compute dose-response curve (DRC) parameters
1. Only attempt on cell line x drug combinations with more than 4 samples (doses)
2. Fit curve to dose and LFC (repeat with LFC.cb)
3. **Resulting data:** `DRC_TABLE`, `DRC_TABLE_cb`

### 8) Collapse LFC table across replicates
1. Take median LFC for each condition
2. **Resulting data:** `LFC_COLLAPSED_TABLE`

### 9) Generate DRC plots for each compound
1. Plot dose-response curves for each cell line x drug
2. **Resulting data:** creates a PDF for each drug in the set containing DRC for each cell line

---

To run biomarker analyses, use `MTS_Analysis.R`, passing a folder with the data output from this pipeline. Output and tests done in the biomarker analysis portion are described in `analysis_info.md`
