# Run after audit-run.R; all subprocess output stays in the supplied directory.
# Rscript tools/test-cli.R CANDIDATE_LIBRARY REFERENCE.rds OUTPUT_DIRECTORY
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
.libPaths(c(normalizePath(args[1]), .libPaths()))
library(phylter)
data(carnivora)
reference <- readRDS(args[2])$results$default
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[3])
exe <- system.file("exec", "phylter", package = "phylter")
stopifnot(nzchar(exe))
# Child R processes use exactly the same library search path.
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
run <- function(options, status = 0L, label) {
  stdout <- file.path(out, paste0(label, ".stdout"))
  stderr <- file.path(out, paste0(label, ".stderr"))
  got <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(exe, options)), stdout = stdout, stderr = stderr))
  if (got != status) stop(label, ": unexpected exit ", got, "\n", paste(readLines(stderr), collapse = "\n"))
  invisible(list(stdout = readLines(stdout), stderr = readLines(stderr)))
}
stopifnot(any(grepl("Usage:", run("--help", label = "help")$stdout)))
run("--version", label = "version")
input <- file.path(out, "carnivora trees.nwk")
ape::write.tree(carnivora, file = input, digits = 17)
prefix <- file.path(out, "analysis result")
logs <- run(c("--trees", input, "--out", prefix, "--quiet", "--save-rds", "--report"), label = "full")
stopifnot(!length(logs$stdout),
          !any(grepl("Initial score|New score|STOPPING OPTIMIZATION", logs$stderr)))
new <- readRDS(paste0(prefix, ".rds"))
# ape's file naming conventions may differ from names on the original R list.
parsed <- ape::read.tree(input)
ids <- names(parsed)
if (is.null(ids)) ids <- as.character(seq_along(parsed))
expected <- reference$Final$Outliers
expected[, 1] <- ids[match(expected[, 1], names(carnivora))]
stopifnot(identical(unname(new$Final$Outliers), unname(expected)))
stopifnot(isTRUE(all.equal(new$Final$AllOptiScores, reference$Final$AllOptiScores, tolerance = 1e-8)))
tab <- read.delim(paste0(prefix, ".outliers.tsv"), colClasses = "character")
stopifnot(identical(names(tab), c("gene", "species")),
          identical(unname(as.matrix(tab)), unname(expected)))
stopifnot(nrow(read.delim(paste0(prefix, ".discarded.tsv"))) == 0L,
          file.info(paste0(prefix, ".pdf"))$size > 0)
before <- tools::md5sum(paste0(prefix, ".outliers.tsv"))
run(c("--trees", input, "--out", prefix), status = 1L, label = "refuse_overwrite")
stopifnot(identical(before, tools::md5sum(paste0(prefix, ".outliers.tsv"))))
run(c("--trees", input, "--out", file.path(out, "bad"), "--k", "oops"), 1L, "invalid_number")
run(c("--unknown"), 1L, "unknown_flag")
run(c("--trees"), 1L, "missing_value")
single <- file.path(out, "single.nwk")
ape::write.tree(carnivora[[1]], file = single)
run(c("--trees", single, "--out", file.path(out, "single")), 1L, "single_tree")
directory <- file.path(out, "directory input")
dir.create(directory, showWarnings = FALSE)
for (i in seq_along(carnivora)) ape::write.tree(carnivora[[i]], digits = 17,
  file = file.path(directory, sprintf("gene%03d.treefile", i)))
dir.prefix <- file.path(out, "directory-result")
run(c("--trees", directory, "--out", dir.prefix, "--quiet"), label = "directory")
tab <- read.delim(paste0(dir.prefix, ".outliers.tsv"), colClasses = "character")
expected <- reference$Final$Outliers
expected[, 1] <- sprintf("gene%03d.treefile", match(expected[, 1], names(carnivora)))
stopifnot(identical(unname(as.matrix(tab)), unname(expected)))
cat("PASS CLI: reference outliers/scores, PDF, TSVs, paths with spaces, directory IDs, errors, overwrite protection\n")
