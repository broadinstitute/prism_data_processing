# packages required for pipeline
options(repos=structure(c(CRAN="http://cran.r-project.org")))

install.packages(c("tidyverse",
                   "useful",
                   "magrittr",
                   "scam",
                   "dr4pl",
                   "readr",
                   "ranger",
                   "sandwich",
                   "BiocManager",
                   "ashr",
                   "corpcor",
                   "lmtest",
                   "PRROC")
                 )

BiocManager::install()
BiocManager::install(c("limma", "sva", "WGCNA"))
