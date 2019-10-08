# FAQs
Frequently asked questions about MTS data processing and reports. This document is continuously updated.

### Is the data used to generate the report public?
All data used for biomarker analysis is available at [depmap.org](https://www.depmap.org). While the exact data is only available to internal collaborators, equivalent public datasets can be found there.

### Why are lineage and mutations analyzed differently than other biomarkers?
A cell line being of a certain lineage or having a mutation is a discrete variable (for example, it either has a given mutation or does not). Because of this, these variables should not tested by correlation with a given outcome. Instead we do a t-test comparing the mean outcome for the set having that feature and the overall mean to see if they are different. Effect size denotes the size of that difference, while q denotes significance. This is analogous to correlation and q for continuous variables (for example, RNA expression).

### How do I know the directionality of the relationship between a given feature and profile in the multivariate analysis?

For the elastic net models, the directionality can be determined by the sign of the mean coefficient, with negative values corresponding to decreased profile values (for example, lower AUC and therefore increased sensitivity). Unfortunately it is not possible to extract this information from the random forest models, but you can typically guess the relationship based on univariate analysis.
