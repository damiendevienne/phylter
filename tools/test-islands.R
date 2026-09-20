# Exhaustively compare the nested old/new island helpers without numerical noise.
# Rscript tools/test-islands.R REFERENCE_SOURCE_DIRECTORY
args <- commandArgs(TRUE)
helper <- function(path) {
  env <- new.env()
  sys.source(path, env)
  for (expr in as.list(body(env$detect.outliers))[-1]) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("detect.island"))) {
      eval(expr, env)
      return(env$detect.island)
    }
  }
  stop("Cannot find island helper")
}
old <- helper(file.path(args[1], "R", "detect.outliers.R"))
new <- helper("R/detect.outliers.R")
for (n in 1:12) {
  for (mask in 0:(2^n - 1)) {
    x <- setNames(as.integer(intToBits(mask))[seq_len(n)], paste0("s", seq_len(n)))
    stopifnot(identical(old(x), new(x)))
  }
}
set.seed(91)
for (i in 1:100) {
  x <- setNames(sample(0:1, 50, TRUE), paste0("s", 1:50))
  stopifnot(identical(old(x), new(x)))
  names(x)[1] <- "out"
  stopifnot(identical(old(x), new(x)))
}
cat("PASS 8,190 exhaustive island patterns plus 200 randomized/sentinel patterns\n")
