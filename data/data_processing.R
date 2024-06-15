# loading in necessary libraries
library(devtools)
library(SummarizedExperiment) # for working with SE objects
library('bigPint') # to get dataset
library(dplyr)
library(tidyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggthemes)


# original data source: https://rdrr.io/bioc/bigPint/man/se_soybean_cn_sub.html
#set working directory
#setwd('#######')


# adding metadata
se_soybean_cn_sub$sample <- colnames(se_soybean_cn_sub)
se_soybean_cn_sub$growth_stage <- c(rep("early",3),rep("middle",3),rep("late",3))
se_soybean_cn_sub$replicate <- rep(c(1,2,3),3)

se_soybean_cn_sub$sample <- colnames(se_soybean_cn_sub)

# writing metadata csv
write.csv(data.frame(metadata(se_soybean_cn_sub)),"metadata_input.csv", row.names = FALSE)

# writing actual data csv for normalized counts data
write.csv(data.frame(assay(se_soybean_cn_sub)),"normalizedcountsmatrix_input.csv", row.names = TRUE)

#writing differential expression data
write.csv(data.frame(rowData(se_soybean_cn_sub)),"differentialexpression_input.csv", row.names = FALSE)

