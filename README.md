## Lilikoi is a novel tool for personalized pathway analysis of metabolomics data.

Lilikoi computes the pathway deregulation score for a given set of metabolites, selects the pathways with the highest mutual information and then uses them to build a classifier.

"Lilikoi: an R package for personalized pathway-based classification modeling using metabolomics data. F. Alakwaa, S. Huang, and L. Garmire (2018) \doi{10.1101/283408}."

## Installation

```
install.packages("lilikoi")

# Or for the latest dev version:
devtools::install_github("lanagarmire/lilikoi")
```

## Example

```
# library(lilikoi)

filename <- system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi")
metaboliteMeasurements <- read.csv(file = filename, check.names = FALSE, row.names = 1)
metaboliteNames <- colnames(metaboliteMeasurements)[-1]
clinicalFactorsData <- read.csv(file = system.file("extdata", "plasma_breast_cancer_Meta.csv",
  package = "lilikoi"))

# The below lines shrink the dataset for faster test runs. Remove them to operate on
# full dataset
metaboliteMeasurements <- metaboliteMeasurements[, 1:20]
metaboliteNames <- colnames(metaboliteMeasurements)[-1]

metabolitePathwayTable <- lilikoi.metab_to_pathway(metaboliteNames, "name")

# We use a subset of the database to speed up tests.
# Swap the comments on the below two lines to run on the full database.
# PDSmatrix <- lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable)
PDSmatrix <- lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable,
  lilikoi::data.smpdb[1:25,])


significantPathways <- lilikoi.select_pathways(PDSmatrix, metaboliteMeasurements,
  threshold = 0.42, method = "gain")

mlResults <- lilikoi.machine_learning(PDSmatrix, metaboliteMeasurements$Label,
  significantPathways)

finalModel <- lilikoi.adjust_model(mlResults$mlResults, PDSmatrix, significantPathways,
  metaboliteMeasurements, clinicalFactorsData, factors = c("Age", "Race"))
```

## Updating the External Databases

Lilikoi depends on data from HMDB, SMPDB, and MetaboAnalyst. This library ships with the latest data as
of the date of publication. To update to the latest data from these sources, load and run the
`lilikoi.update_database()` method found in the `lilikoi.update_database.r` file.

Warning: the datasets are large (>5GB) and this step may take greater than 20 minutes.

# Built By

*   Fadhl Alakwaa https://github.com/FADHLyemen
*   Sijia Huang https://github.com/scarlettcanny

# More Examples

*   Shiny Version: http://lilikoi.garmiregroup.org
*   https://github.com/lanagarmire/lilikoi/blob/master/lilikoi_example.ipynb
*   https://github.com/lanagarmire/lilikoi/blob/master/lilikoiExample.r
*   https://mybinder.org/v2/gh/FADHLyemen/lilikoi_Fadhl/master
