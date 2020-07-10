# TIME
**Description:** Use gene expression data to get the corresponding TIME phenotypes and TIMEscore

**Note:** Make sure that the data you enter is a matrix, and rows are samples and columns are genes.

## Download
library(devtools)

install_github("Zaoqu-Liu/TIME")

## An example is as follows
if(!require("GSVA")) BiocManager::install("GSVA",update = F,ask = F)

if(!require("clusterRepro")) install.packages("clusterRepro",update = F,ask = F)

library(TIME)

data('data')

result=time_phenotype(data2)
