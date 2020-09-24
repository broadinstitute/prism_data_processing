# Column Headers

## Data

### logMFI and logMFI\_NORMALIZED
- **pert\_mfc\_id:** the Broad ID for the compound
- **profile\_id:** a concatenation of replicate and well
- **rid:** a unique identifier for this cell line
- **ccle\_name:** the name of the cell line
- **pool\_id:** the PRISM pool the cell line is in
- **culture:** the PRISM cell set the cell line is in
- **prism\_replicate:** the replicate this data point was obtained from
- **pert\_type:** perturbation type, denotes whether this was a positive control, negative control, or treatment condition
- **pert\_dose:** the numeric dose of this condition (µM)
- **pert\_idose:** a string version of pert_dose
- **pert\_name:** the name of the compound
- **pert\_well:** the plate well the cell line was in
- **pert\_time:** the length of the assay (in hours)
- **logMFI:** the log2 mean-fluorescent intensity (Luminex readout)
- **LMFI:** normalized logMFI

### LFC\_TABLE
- **pert\_mfc\_id:** the Broad ID for the compound
- **pert\_name:** the name of the compound
- **prism\_replicate:** the replicate this data point was obtained from
- **culture:** the PRISM cell set the cell line is in
- **rid:** a unique identifier for this cell line
- **LFC:** the log-fold change of this condition versus DMSO
- **pert\_type:** perturbation type, denotes whether this was a positive control, negative control, or treatment condition
- **ccle\_name:** the name of the cell line
- **pert\_dose:** the numeric dose of this condition (µM)
- **pert\_well:** the plate well the cell line was in
- **pert\_time:** the length of the assay (in hours)
- **pool\_id:** the PRISM pool the cell line is in
- **profile\_id:** a concatenation of replicate and well
- **pert\_idose:** a string version of pert_dose
- **compound\_plate:** the plate of the sample
- **LFC.cb:** the COMBAT corrected log-fold change (this is used for the dose-response curves, see the report for more information)

### LFC\_COLLAPSED\_TABLE
- **ccle\_name:** the name of the cell line
- **culture:** the PRISM cell set the cell line is in
- **pert\_name:** the name of the compound
- **pert\_mfc\_id:** the Broad ID for the compound
- **pert\_dose:** the numeric dose of this condition (µM)
- **pert\_idose:** a string version of pert_dose
- **compound\_plate:** the plate of the sample
- **pert\_time:** the length of the assay (in hours)
- **LFC:** the log-fold change of this condition versus DMSO
- **LFC.cb:** the COMBAT corrected log-fold change (this is used for the dose-response curves, see the report for more information)

### DRC\_TABLE
- **min\_dose, max\_dose:** the minimum and maximum doses at which a compound was tested
- **upper\_limit, lower\_limit:** the asymptotic maximum and minimum of the curve
- **ec50:** the dose at which the curve reaches the value between its upper and lower limit
- **slope:** the slope of the curve at the EC50
- **auc:** area under the curve
- **log2.ic50:** log2 of the IC50 value (the point at which the cell line reaches 50% viability)
- **mse:** the mean-squared error of the curve
- **R2:** the R-squared of the curve
- **ccle\_name:** the name of the cell line
- **culture:** the PRISM cell set the cell line is in
- **pert\_mfc\_id:** the Broad ID for the compound
- **pert\_name:** the name of the compound
- **pert\_time:** the length of the assay (in hours)


## Results

### SSMD\_TABLE
- **prism\_replicate:** the replicate this data point was obtained from
- **ccle\_name:** the name of the cell line
- **pool\_id:** the PRISM pool the cell line is in
- **culture:** the PRISM cell set the cell line is in
- **ctl\_vehicle\_md:** median logMFI value for negative controls
- **trt\_poscon\_md:** median logMFI value for positive controls
- **ctl\_vehicle\_mad:** median absolute deviation of negative controls
- **trt\_poscon\_mad:** median absolute deviation of positive controls
- **ssmd:** strictly-standardizedmean difference of positive and negative controls
- **nnmd:** null-normalized mean difference of positive and negative controls
- **rid:** a unique identifier for this cell line
- **error\_rate:** the error rate of a perfect threshold classifier distinguishing positive and negative controls
- **pass:** whether this cell line passed QC or not
- **pert\_time:** the length of the assay (in hours)
- **compound\_plate:** the plate of the sample
- **n.rep:** number of passing replicates for this cell line

### continuous\_associations
- **feature:** the feature that was correlated
- **PosteriorMean:** the adaptive shrinkage moderated effect size estimates
- **PosteriorSD:** the standard deviation of the PosteriorMean
- **qvalue:** the false-discovery rate of the PosteriorMean
- **coef:** the correlation coefficient between the feature and the response variable
- **q.val:** the false-discovery rate of the coef
- **rank:** the rank of the strength of this correlation
- **pert\_mfc\_id:** the Broad ID for the compound
- **pert\_name:** the name of the compound
- **dose:** the response variable (either log-fold change at a particular dose or AUC or IC50)
- **feature\_type:** the type of feature
- **pert\_time:** the length of the assay (in hours)

### discrete\_associations
- **feature:** the feature that was correlated
- **effect\_size:** the difference in mean response of cell lines with this feature versus all others
- **t:** the t-test result
- **p:** the p-value of the t-test
- **q:** the false-discovery rate corrected p-value
- **pert\_mfc\_id:** the Broad ID for the compound
- **pert\_name:** the name of the compound
- **dose:** the response variable (either log-fold change at a particular dose or AUC or IC50)
- **feature\_type:** the type of feature
- **pert\_time:** the length of the assay (in hours)

### Model\_table
- **MSE:** the mean-squared error of the model
- **MSE.se:** the standard error of the mean-squared error
- **R2:** the R-squared of the model
- **PearsonScore:** the Pearson score of the model
- **type:** the type of model fit
- **pert\_mfc\_id:** the Broad ID for the compound
- **pert\_name:** the name of the compound
- **pert\_idose:** the response variable (either log-fold change at a particular dose or AUC or IC50)
- **model:** the dataset used for the model
- **pert\_time:** the length of the assay (in hours)

### RF\_table
- **feature:** the omic feature being considered
- **RF.imp.mean:** the mean importance of that feature across cross-validation steps
- **RF.imp.sd:** the standard deviation of the importance
- **RF.imp.stability:** the fraction of models using that feature
- **pert\_mfc_id:** the Broad ID for the compound
- **pert\_name:** the name of the compound
- **pert\_idose:** the response variable (either log-fold change at a particular dose or AUC or IC50)
rank
- **model:** the dataset used for the model
- **pert\_time:** the length of the assay (in hours)
