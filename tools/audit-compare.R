# Rscript tools/audit-compare.R REFERENCE.rds CANDIDATE.rds
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
old <- readRDS(args[1])
new <- readRDS(args[2])
near <- function(a, b, label) {
  result <- all.equal(a, b, tolerance = 1e-8)
  if (!isTRUE(result)) stop(label, ": ", paste(result, collapse = "; "))
  # all.equal uses an aggregate scale: also bound every finite numeric entry.
  check_entries <- function(x, y) {
    if (is.list(x)) {
      for (i in seq_along(x)) check_entries(x[[i]], y[[i]])
    } else if (is.numeric(x)) {
      finite <- is.finite(x) & is.finite(y)
      if (any(abs(x[finite] - y[finite]) > 1e-9 + 1e-8 * pmax(abs(x[finite]), abs(y[finite]))))
        stop(label, ": elementwise numerical tolerance exceeded")
    }
  }
  check_entries(a, b)
}
projection <- function(a, b, label) {
  stopifnot(identical(dim(a$F), dim(b$F)), identical(dimnames(a$F), dimnames(b$F)))
  # Eigenvector signs (and bases within tied eigenspaces) are arbitrary.
  near(tcrossprod(a$F), tcrossprod(b$F), paste(label, "F geometry"))
  stopifnot(identical(names(a$PartialF), names(b$PartialF)))
  for (i in seq_along(a$PartialF)) {
    near(tcrossprod(a$PartialF[[i]], a$F), tcrossprod(b$PartialF[[i]], b$F),
         paste(label, "partial/compromise geometry", i))
    near(tcrossprod(a$PartialF[[i]]), tcrossprod(b$PartialF[[i]]),
         paste(label, "partial geometry", i))
  }
}
stage <- function(a, b, label) {
  for (field in c("WR", "RV", "weights", "compromise", "matrices", "mat.data"))
    near(a[[field]], b[[field]], paste(label, field))
  projection(a, b, label)
}
for (name in setdiff(names(old$results), "kernel")) {
  a <- old$results[[name]]
  b <- new$results[[name]]
  stopifnot(identical(class(a), class(b)))
  if (name == "initial") stage(a, b, name) else {
    stage(a$Initial, b$Initial, paste(name, "initial"))
    stage(a$Final, b$Final, paste(name, "final"))
    for (field in c("Outliers", "CELLSREMOVED", "CompleteOutliers", "Discarded", "species.order"))
      if (!identical(a$Final[[field]], b$Final[[field]])) stop(name, ": different ", field)
    stopifnot(identical(a$DiscardedGenes, b$DiscardedGenes))
    near(a$Final$AllOptiScores, b$Final$AllOptiScores, paste(name, "optimization trajectory"))
  }
  cat("PASS", name, "\n")
}
a <- old$results$kernel
b <- new$results$kernel
near(a$imputed, b$imputed, "synthetic imputation")
near(a$wr, b$wr, "synthetic WR")
near(a$medcouple, b$medcouple, "medcouple edge cases")
for (field in c("alpha", "lambda", "RVmat", "compromise", "quality", "matrices.dblcent"))
  near(a$distatis[[field]], b$distatis[[field]], paste("synthetic", field))
projection(a$distatis, b$distatis, "synthetic")
cat("PASS synthetic kernels\n")
for (name in names(old$timings)) {
  x <- median(old$timings[[name]])
  y <- median(new$timings[[name]])
  cat(sprintf("%-24s reference %.4fs candidate %.4fs speedup %.2fx\n", name, x, y, x/y))
}
print(rbind(reference_allocated_bytes = old$allocated, candidate_allocated_bytes = new$allocated))
