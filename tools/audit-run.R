# Run from the repository root, in separate R processes for each installation.
# Rscript tools/audit-run.R LIBRARY OUTPUT_PREFIX [REPEATS]
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
.libPaths(c(normalizePath(args[1]), .libPaths()))
library(phylter)
data(carnivora)
repeats <- if (length(args) > 2L) as.integer(args[3]) else 3L
prefix <- args[2]
set.seed(20260920)
missing <- lapply(carnivora, function(tr) ape::drop.tip(tr, sample(tr$tip.label, 5)))
matrices <- trees2matrices(carnivora)
zero <- matrices[[1]] * 0
cases <- list(
  default = list(X = carnivora),
  nodal = list(X = carnivora, distance = "nodal"),
  mean = list(X = carnivora, Norm = "mean"),
  unscaled = list(X = carnivora, Norm = "none"),
  no_islands = list(X = carnivora, test.island = FALSE),
  parallel = list(X = carnivora, parallel = TRUE),
  column_normalization = list(X = carnivora, normalizeby = "col"),
  no_wr_normalization = list(X = carnivora, normalizeby = "none"),
  gene_threshold = list(X = carnivora, k2 = 1.5),
  missing_taxa = list(X = missing),
  discarded = list(X = c(list(zero_gene = zero), matrices)),
  initial = list(X = carnivora, InitialOnly = TRUE)
)
results <- list()
timings <- list()
for (name in names(cases)) {
  gc()
  set.seed(42)
  elapsed <- system.time({
    results[[name]] <- do.call(phylter, modifyList(list(parallel = FALSE, verbose = FALSE), cases[[name]]))
  })[["elapsed"]]
  timings[[name]] <- elapsed
  cat(name, elapsed, "seconds\n")
}
for (i in seq_len(repeats - 1L)) {
  gc()
  set.seed(42)
  timings$default <- c(timings$default, system.time(
    invisible(capture.output(phylter(carnivora, parallel = FALSE, verbose = FALSE)))
  )[["elapsed"]])
}
# Moderate synthetic kernels: independently vary the species count beyond 53.
set.seed(123)
large <- lapply(seq_len(80), function(i) {
  m <- as.matrix(dist(matrix(rnorm(200 * 6), 200)))
  dimnames(m) <- list(paste0("s", seq_len(200)), paste0("s", seq_len(200)))
  m
})
names(large) <- paste0("g", seq_along(large))
incomplete <- lapply(large, function(m) {
  keep <- sample(nrow(m), 170)
  m[keep, keep]
})
bench <- function(fun) {
  vapply(seq_len(repeats), function(i) { gc(); system.time(fun())[["elapsed"]] }, numeric(1))
}
timings$distatis_200x80 <- bench(function() DistatisFast(large, parallel = FALSE))
timings$imputation_200x80 <- bench(function() impMean(incomplete))
ds <- DistatisFast(large, parallel = FALSE)
timings$wr_200x80 <- bench(function() for (i in 1:20) Dist2WR(ds))
results$kernel <- list(distatis = ds, imputed = impMean(incomplete), wr = Dist2WR(ds))
mc.inputs <- list(single = 1, pair = c(1, 2), constant = rep(1, 200),
                  tied = c(rep(0, 60), 1:40), symmetric = -50:50,
                  skewed = exp(rnorm(201)), short = rnorm(20))
results$kernel$medcouple <- lapply(mc.inputs, function(x)
  lapply(list(NULL, FALSE, TRUE), function(reflect) medcouple(x, reflect)))
# Report allocations separately: this is cumulative R allocation, not peak RSS.
allocation <- function(fun, suffix) {
  path <- paste0(prefix, "-", suffix, ".Rprofmem")
  gc()
  Rprofmem(path)
  tryCatch(fun(), finally = Rprofmem(NULL))
  lines <- readLines(path)
  sum(suppressWarnings(as.numeric(sub(" .*", "", lines))), na.rm = TRUE)
}
allocated <- c(
  distatis = allocation(function() DistatisFast(large, parallel = FALSE), "distatis"),
  imputation = allocation(function() impMean(incomplete), "imputation")
)
Rprof(paste0(prefix, ".Rprof"), interval = 0.01)
invisible(capture.output(phylter(carnivora, parallel = FALSE, verbose = FALSE)))
Rprof(NULL)
saveRDS(list(results = results, timings = timings, allocated = allocated,
             session = sessionInfo(), profile = summaryRprof(paste0(prefix, ".Rprof"))),
        paste0(prefix, ".rds"))
print(timings)
print(allocated)
print(head(summaryRprof(paste0(prefix, ".Rprof"))$by.total, 15))
