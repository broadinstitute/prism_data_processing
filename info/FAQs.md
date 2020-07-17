# FAQs
Frequently asked questions about MTS data processing and reports. This document is continuously updated.

### Which datasets were used for biomarker analysis?

For external collaborators, the latest public datasets available on [depmap.org](https://depmap.org/portal) are used. For internal collaborators, the latest internal datasets are used, they are available on internal versions of the portal.

### Why use error rate instead of SSMD as QC?

We find that error rate is better at determining which cell lines are noisiest. If you wish to use SSMD to QC your data (for example to reflect previous PRISM results) the information is tabulated in `SSMD_TABLE.csv`. You can also find other metrics like NNMD and normalized logMFI values for positive and negative controls.

### Why are mutation and lineage separate from other biomarkers?

Lineage and mutation are discrete variables that describe a cell line and therefore require different statistical tests than continuous variables to quantify their effect (t-test as opposed to correlation). For more information on how biomarker analyses are done see our [analysis info document](./analysis_info.pdf).

### How can I filter cell lines down to have more "trustworthy" lines?

We only use error rate to determine whether cell lines pass QC (see the [pipeline description](./MTS_pipeline.md) for details). We do this to be as light-handed with the data as possible. To filter further we recommend a few different strategies:

1. Filter by **SSMD**. Typically a cutoff of 2 is used (higher is better).
2. Look at **replicate correlation**. Calculate the correlation between each replicate for each cell line (there are 3 replicates for each) in `LFC_TABLE.csv` and use only cell lines with reasonably good correlations across replicates. It is important to not that cell lines that are killed more tend to have higher replicate correlations so filtering in this way can remove more insensitive than sensitive cell lines.
3. Look at the **$R^2$ and MSE of the dose-response curves** in `DRC_TABLE.csv` and filter for cell lines that have a good curve fit.
4. Calculate the **monotonicity** of each cell line's dose-response curve. Filter out cell lines that do not respond monotonically (or nearly monotonically) to your compound.
