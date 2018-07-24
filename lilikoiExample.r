library(lilikoi)

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
