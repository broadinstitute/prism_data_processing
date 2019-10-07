# FAQs
Frequently asked questions about MTS data processing and reports.

### How do I know the directionality of the relationship between a given feature and profile in the multivariate analysis?

For the elastic net models, the directionality can be determined by the sign of the mean coefficient, with negative values corresponding to decreased profile values (for example, lower AUC and therefore increased sensitivity). Unfortunately it is not possible to extract this information from the random forest models, but you can typically guess the relationship based on univariate analysis.

### Is the data used to generate the report public?
