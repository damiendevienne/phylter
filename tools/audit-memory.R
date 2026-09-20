# Run separately under /usr/bin/time -v: peak RSS includes R and input storage.
# Rscript tools/audit-memory.R LIBRARY distatis|imputation [SPECIES] [GENES]
args <- commandArgs(TRUE)
.libPaths(c(normalizePath(args[1]), .libPaths()))
library(phylter)
n <- if (length(args) >= 3L) as.integer(args[3]) else 400L
g <- if (length(args) >= 4L) as.integer(args[4]) else 150L
set.seed(123)
m <- lapply(seq_len(g), function(i) {
  x <- as.matrix(dist(matrix(rnorm(n * 6), n)))
  dimnames(x) <- list(paste0("s", seq_len(n)), paste0("s", seq_len(n)))
  if (args[2] == "imputation") {
    keep <- sample(n, floor(n * 0.85))
    x <- x[keep, keep]
  }
  x
})
names(m) <- paste0("g", seq_len(g))
gc()
if (args[2] == "distatis") invisible(DistatisFast(m, parallel = FALSE)) else
  if (args[2] == "imputation") invisible(impMean(m)) else stop("Unknown workload")
