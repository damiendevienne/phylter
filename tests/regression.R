library(phylter)
data(carnivora)
fixture <- file.path("fixtures", "carnivora-reference.rds")
if (!file.exists(fixture)) fixture <- file.path("tests", fixture)
reference <- readRDS(fixture)
set.seed(42)
result <- phylter(carnivora, parallel = FALSE, verbose = FALSE)
for (stage in c("Initial", "Final")) {
  expected <- reference[[stage]]
  actual <- result[[stage]]
  for (field in c("WR", "RV", "weights", "compromise"))
    stopifnot(isTRUE(all.equal(expected[[field]], actual[[field]], tolerance = 1e-8)))
  stopifnot(isTRUE(all.equal(expected$FGram, tcrossprod(actual$F), tolerance = 1e-8)))
}
for (field in c("Outliers", "CELLSREMOVED", "CompleteOutliers", "Discarded", "species.order"))
  stopifnot(identical(reference$Final[[field]], result$Final[[field]]))
stopifnot(identical(reference$DiscardedGenes, result$DiscardedGenes),
          isTRUE(all.equal(reference$Final$AllOptiScores, result$Final$AllOptiScores, tolerance = 1e-8)))
