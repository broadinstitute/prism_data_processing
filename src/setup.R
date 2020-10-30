# packages required for pipeline
options(repos=structure(c(CRAN="http://cran.r-project.org")))

install.packages("tidyverse")
install.packages("hdf5r")
install.packages("reshape2")
install.packages("readr")
install.packages("magrittr")
install.packages("dr4pl")
install.packages("data.table")
install.packages("scam")
install.packages("PRROC")
install.packages("gmodels")
install.packages("devtools")
install.packages("BiocManager") # to install limma

# BiocManager packages
BiocManager::install()
BiocManager::install("sva")
BiocManager::install("limma", ask = FALSE)

# CDS github packages
devtools::install_github("broadinstitute/taigr", dependencies=TRUE, force=TRUE)
devtools::install_github("broadinstitute/cdsr_models", dependencies=TRUE, force=TRUE)
