context("test-all.R")

test_that("full pipeline works", {

  filename <- system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi")

  metaboliteMeasurements <- read.csv(file = filename, check.names = F, row.names = 1)
  metaboliteMeasurements <- metaboliteMeasurements[, 1:20]
  metaboliteNames <- colnames(metaboliteMeasurements)[-1]

  metabolitePathwayTable <- lilikoi.metab_to_pathway(metaboliteNames, "name")
  PDSmatrix <- lilikoi.get_pd_scores(metaboliteMeasurements, metabolitePathwayTable, lilikoi::data.smpdb[1:23,], maxit = 1)

  significantPathways <- lilikoi.select_pathways(PDSmatrix, metaboliteMeasurements, threshold = 0.42, method = "gain")

  #mlResults <- lilikoi.machine_learning(PDSmatrix, metaboliteMeasurements$Label, significantPathways)
  #result <- mlResults$res

  #clinicalFactorsData <- read.csv(file = system.file("extdata", "plasma_breast_cancer_Meta.csv", package = "lilikoi"))
  #model = lilikoi.adjust_model(result, PDSmatrix, significantPathways, metaboliteMeasurements, clinicalFactorsData, factors = c("Age", "Race"))

  # Todo: add more checks above
  expect_equal(1, 1)
})
