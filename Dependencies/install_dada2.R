# If install fails, try changing "BiocVersion" variable in line 13 based on your R version:
# R version = BiocManager Version
# 4.0.2+ = 3.12
# 4.0 = 3.11
# 3.6 = 3.10
# 3.5 = 3.8
# 3.4 = 3.6
# 3.3 = 3.4
# 3.2 = 3.2

if (R.version$major == "4"){
  if (R.version$minor == "1.0"){
    BiocVersion <- 3.13
  } else if (R.version$minor %in% c("0.5","0.4","0.3","0.2","0.1","0.0")){
    BiocVersion <- 3.12
  }
} else if (R.version$major == "3"){
  if (R.version$minor %in% c("6.3","6.2","6.1","6.0")){
    BiocVersion <- 3.10
  } else if (R.version$minor %in% c("5.3","5.2","5.1","5.0")){
    BiocVersion <- 3.8
  }
}

### Install ###
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", quiet = TRUE)
BiocManager::install(version = BiocVersion, ask = FALSE)
if (!require("dada2", quietly = TRUE)) BiocManager::install("dada2", version = "3.12", ask = FALSE)
if (!require("BioStrings", quietly = TRUE)) BiocManager::install("Biostrings", ask = FALSE)
if (!require("ShortRead", quietly = TRUE)) BiocManager::install("ShortRead", ask = FALSE)
if (!require("gridExtra", quietly = TRUE)) install.packages("gridExtra", quiet = TRUE)
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2", quiet = TRUE)
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2", quiet = TRUE)
if (!require("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", quiet = TRUE)
