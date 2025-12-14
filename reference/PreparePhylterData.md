# Prepare data for phylter analysis

Prepare datasets for the `phylter` function. Detects possible issues,
discards genes if necessary, imputes missing data if any, and reorders
row- and col-names. For internal usage mostly.

## Usage

``` r
PreparePhylterData(
  X,
  bvalue = 0,
  distance = "patristic",
  Norm = "median",
  Norm.cutoff = 0.001,
  gene.names = NULL,
  verbose = TRUE
)
```

## Arguments

- X:

  A list of phylogenetic trees (phylo object) or a list of distance
  matrices. Trees can have different number of leaves and matrices can
  have different dimensions. If this is the case, missing values are
  imputed.

- bvalue:

  If X is a list of trees, nodes with a support below 'bvalue' will be
  collapsed prior to the outlier detection.

- distance:

  If X is a list of trees, type of distance used to compute the pairwise
  matrices for each tree. Can be "patristic" (sum of branch lengths
  separating tips, the default) or nodal (number of nodes separating
  tips).

- Norm:

  Should the matrices be normalized and how. If "median", matrices are
  divided by their median, if "mean" they are divided by their mean, if
  "none", no normalization if performed. Normalizing ensures that
  fast-evolving (and slow-evolving) genes are not treated as outliers.
  Normalization by median is less sensitive to outlier values but can
  lead to errors if some matrices have a median value of 0. are not
  considered outliers.

- Norm.cutoff:

  Value of the median (if Norm="median") or the mean (if Norm="mean")
  below which matrices are simply discarded from the analysis. This
  prevents dividing by 0, and getting rid of genes that contain mostly
  branches of length 0 and are therefore uninformative anyway.

- gene.names:

  List of gene names used to rename elements in X. If NULL (the
  default), elements are named 1,2,..,length(X).

- verbose:

  If TRUE (the default), messages are written during the filtering
  process to get information of what is happening

## Value

A list of class 'phylter' with the 'Initial' (before filtering) and
'Final' (after filtering) states, or a list of class 'phylterinitial'
only, if InitialOnly=TRUE.

## Examples

``` r
data(carnivora)
# transform trees to a named list of matrices with same dimensions
# and identical row and column orders and names
carnivora.clean<-PreparePhylterData(carnivora)
#> 
#> Number of Genes:    125
#> Number of Species:  53
#> --------

 
```
