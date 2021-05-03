# If install fails, try changing BiocManager version in lines 10-11 based on your R version:
# R version = BiocManager Version
# 4.0.2+ = 3.12
# 4.0 = 3.11
# 3.6 = 3.10
# 3.5 = 3.8
# 3.4 = 3.6
# 3.3 = 3.4
# 3.2 = 3.2

### Install ###
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
BiocManager::install(version = '3.12', ask = FALSE)
if (!require("dada2", quietly = TRUE)) BiocManager::install("dada2", version = "3.12", ask = FALSE)
if (!require("Biostrings", quietly = TRUE)) BiocManager::install("Biostrings", ask = FALSE)
if (!require("ShortRead", quietly = TRUE)) BiocManager::install("ShortRead", ask = FALSE)
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2", quiet = TRUE)
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2", quiet = TRUE)
if (!require("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", quiet = TRUE)
