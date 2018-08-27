## Lilikoi is a novel tool for personalized pathway analysis of metabolomics data. 

## Prerequisites

To install all the required packages without overwriting your installed packages, you can run the below lines:

```
list.of.packages <- c("ggplot2", "caret","devtools", "dplyr","RWeka","infotheo","pROC","reshape2","corrplot", "stringr", "Hmisc", "Matrix", "randomForest", "glmnet", "gbm", "e1071", "pamr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

source("https://bioconductor.org/biocLite.R")
biocLite("pathifier")
```

## Installation
```
devtools::install_github("lanagarmire/lilikoi")
```

If you have a problem installing Rweka as it requires Java, you can reconfigure R from the command line by running the below line:

```
R CMD javareconf
```

# Built By
* Fadhl Alakwaa https://github.com/FADHLyemen
* Sijia Huang  https://github.com/scarlettcanny

# Example Code
https://github.com/lanagarmire/lilikoi/blob/master/lilikoi_example.ipynb
https://mybinder.org/v2/gh/FADHLyemen/lilikoi_Fadhl/master

# To contribute
library(devtools)
