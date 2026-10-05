library(phylter)
data(carnivora)
near <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-10)))

# Dist2WR must preserve names and shape, including one-factor projections.
set.seed(17)
for (nf in c(1L, 2L, 7L)) {
  F <- matrix(rnorm(11 * nf), 11, dimnames = list(paste0("s", 1:11), NULL))
  partial <- list(g1 = F + 0.2, g2 = F * 2, g3 = F - 1)
  expected <- do.call(cbind, lapply(partial, function(x)
    apply((x - F)^2, 1, function(y) sqrt(sum(y)))))
  near(Dist2WR(list(F = F, PartialF = partial)), expected)
}

# A separate centering/weighted-mean identity checks the native accumulator.
m <- PreparePhylterData(carnivora[1:12], verbose = FALSE)$matrices
before <- serialize(m, NULL)
ds <- DistatisFast(m, factorskept = 3, parallel = FALSE)
S <- Reduce(`+`, Map(function(x, w) x * w, ds$matrices.dblcent, ds$alpha))
expected <- outer(diag(S), diag(S), `+`) - 2 * S
dimnames(expected) <- dimnames(m[[1]])
near(ds$compromise, expected)
stopifnot(identical(before, serialize(m, NULL)))

# Missing taxa and reordered labels: compare with a direct 3D array mean.
missing <- m
missing[[1]] <- m[[1]][-c(1, 3), -c(1, 3)]
missing[[2]] <- m[[2]][53:1, 53:1]
species <- unique(unlist(lapply(missing, rownames)))
expanded <- lapply(missing, function(x) {
  y <- matrix(NA_real_, length(species), length(species), dimnames = list(species, species))
  y[rownames(x), colnames(x)] <- x
  y
})
arr <- array(unlist(expanded), c(length(species), length(species), length(missing)))
average <- apply(arr, c(1, 2), mean, na.rm = TRUE)
expected <- lapply(expanded, function(x) { x[is.na(x)] <- average[is.na(x)]; x })
near(impMean(missing), expected)

# Exercise the native early-return paths that previously leaked workspace.
for (i in 1:100) {
  stopifnot(as.numeric(medcouple(c(1, 2))) == 0)
  stopifnot(is.finite(as.numeric(medcouple(rep(1, 200)))))
}
set.seed(7)
x <- rnorm(201)
near(as.numeric(medcouple(x, TRUE)), -as.numeric(medcouple(-x, TRUE)))
