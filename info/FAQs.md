# FAQs
Frequently asked questions about MTS data processing and reports. This document is continuously updated.

### Which datasets were used for biomarker analysis?

### Why use error rate instead of SSMD as QC?
We find that error rate is better at determining which cell lines are noisiest. If you wish to use SSMD to QC your data (for example to reflect previous PRISM results) the information is tabulated in `SSMD_TABLE.csv`.

### Why are mutation and lineage separate from other biomarkers?
Lineage and mutation are discrete variables that describe a cell line and therefore require different statistical tests than continuous variables to quantify their effect (t-test as opposed to correlation).
