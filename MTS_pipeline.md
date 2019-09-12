# MTS Processing Pipeline

R script to process Level 2 MTS data. This documents describes steps executed by the script.

- **Input** is logMFI values associated with 3 replicates of cell line treatments with drugs and controls.
- **Output** is files describing dose-response of each cell line to each drug, mean-fold change of cell lines in response to drugs, and statistical significance metrics.
- **To execute:** change `data_path` variable to match local path of unprocessed data.

### 1) Load Data
1. Read .csv of logMFI (mean fluorescence increase) data for each sample
2. Split into PR300 and PR500 datasets
3. Create reference tables containing only control barcodes (pool ID is "CTLBC")
4. **Resulting data:** `PR300`, `PR500`, `PR300_barcodes`, `PR500_barcodes`

### 2) Normalize control vehicle samples (negative controls)
1. Calculate median logMFI of each well across replicates
2. Take the median of calculated median to generate overall median
3. Subtract well median from sample value and add back overall median (result is normalized samples with respect to overall median)
4. **Resulting data:** `PR300_normalized` (just negative controls), `PR500_normalized` (just negative controls)

### 3) Normalize all samples using controls
1. Fit a curve to normalized control data (logMFI ~ normalized logMFI)
2. Use curve equation to normalize non-control data
3. Recombine with original data replacing logMFI with normalized values
4. **Resulting data:** `PR300_normalized` (all data), `PR500_normalized` (all data)

### 4) Generate statistics table
1. For each well:
 1. Calculate median difference (MD)
 2. Calculate median absolute deviation (MAD)
 3. Calculate strictly standardized mean difference (SSMD)
 4. Calculate null-normalized mean difference (NNMD)
 5. Generate calculate error rate of logMFI distributions
2. Remove samples with an error rate > 0.05
3. **Resulting data:** `SSMD_TABLE`

### 5) Calculate log-fold changes (LFC)
1. Only use samples with error rate < 0.05
2. Wells must be present or have a low enough error rate on at least 2 replicates
3. $LFC_{sample} = logMFI_{sample} - median(logMFI_{controls})$
4. **Resulting data:** `LFC_TABLE`

### 6) Correct for pool effects using ComBat
1. Generates a corrected LFC for each sample (`LFC.cb`) and stores in new table (separate from uncorrected)
2. **Resulting data:** `LFC_TABLE` (with new LFC.cb value)

### 7) Compute dose-response curce (DRC) parameters
1. Only attempt on cell line x drug combinations with more than 4 samples (doses)
2. Fit curve to dose and LFC (repeat with LFC.cb)
3. **Resulting data:** `DRC_TABLE`, `DRC_TABLE_cb`

### 8) Collapse LFC table across replicates
1. Take median LFC for each condition
2. **Resulting data:** `LFC_COLLAPSED_TABLE`

### 9) Generate DRC plots for each compound
1. Plot dose-response curves for each cell line x drug
2. **Resulting data:** creates a PDF for each drug in the set containing DRC for each cell line
