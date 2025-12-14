# Median normalization of 2D matrix by row or by colomn

This function normalizes the 2WR matrix (or any 2D matrix) according to
the species (rows) or to the genes (columns), by dividing each row or
each column by its median.

## Usage

``` r
normalize(mat, what = "none")
```

## Arguments

- mat:

  A matrix

- what:

  Character string indicating whether the matrix should be normalized
  and how. If what="none", the matrix is not normalized (the default),
  if what="species", the matrix is normalized so that the difference
  between species is increased, and if what="genes", the matrix is
  normalized so that the difference between genes is increased.
  Normalization consists in dividing either each row or each columns by
  its median.

## Value

A normalized matrix

## Examples

``` r
# random matrix
x<-matrix(rnorm(270), nrow=9, ncol=14)
#> Warning: data length [270] is not a sub-multiple or multiple of the number of columns [14]

# normalize by row
x1<-normalize(x, "genes")

# normalize by column
x2<-normalize(x, "species")


```
