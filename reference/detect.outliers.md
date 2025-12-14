# Detection of outliers in 1D and 2D data

Functions to detect outliers, in matrices or in arrays.

## Usage

``` r
detect.outliers(X, k = 3, test.island = TRUE, normalizeby = "row")

detect.outliers.array(arr, nbsp, k = 3)
```

## Arguments

- X:

  2D matrix (gene x species) obtained with the [`Dist2WR()`](Dist2WR.md)
  function.

- k:

  strength of outlier detection. High values (typically 3) correspond to
  less outliers detected than with lower ones (e.g. 1.5).

- test.island:

  should islands of outliers be treated as such. If TRUE (the default),
  only the highest value in an island of outliers is removed. This
  prevents erroneous outliers due to hicthiking effect to be removed.

- normalizeby:

  Should the 2D matrix be normalized prior to outlier detection, and
  how. Can be "row" (the default),"col" or "none". Normalization is done
  by dividing columns or rows by their median.

- arr:

  Array of values, typically the weight of each gene matrix (alpha
  values).

- nbsp:

  Number of species in the analysis

## Value

`detect.outliers`: A matrix with outliers detected in the 2D matric.
Each row `x` contains the gene (`x[1]`) where the species (`x[2]`) is an
outlier.

`detect.outliers.array`: An array listing the outliers detected (if any)

## Details

These functions detect outliers either in matrices or in arrays. For the
method to be adapted to skewed data, as is the case here, the outlier
detection method used is the adjusted Tukey proposed by Hubert and
Vandervieren (2008).

## Functions

- `detect.outliers()`: detect outliers in 2D matrix

- `detect.outliers.array()`: detects outliers in 1D array

## References

de Vienne D.M., Ollier S. et Aguileta G. (2012) Phylo-MCOA: A Fast and
Efficient Method to Detect Outlier Genes and Species in Phylogenomics
Using Multiple Co-inertia Analysis. Molecular Biology and Evolution 29 :
1587-1598.

Hubert, M. and Vandervieren, E. (2008). An adjusted boxplot for skewed
distributions. Computational Statistics and Data Analysis, 52,
5186-5201.

## Examples

``` r
# Get the initial gene x species matrix
# from the carnivora dataset
data(carnivora) 
mat <- phylter(carnivora, InitialOnly = TRUE, parallel = FALSE)$WR
#> 
#> Number of Genes:    125
#> Number of Species:  53
#> --------

# detect outliers in this matrix
outliers<-detect.outliers(mat)
outliers$cells # matrix where each row represents one cell in the input matrix
#>       [,1] [,2]
#>  [1,]    2    6
#>  [2,]    2   33
#>  [3,]    2   42
#>  [4,]   38    4
#>  [5,]   38   53
#>  [6,]   38   30
#>  [7,]   52   48
#>  [8,]   55   28
#>  [9,]   56   31
#> [10,]   60    7
#> [11,]   60   35
#> [12,]   60   38
#> [13,]   60    2
#> [14,]   64   33
#> [15,]   76    3
#> [16,]   76    7
#> [17,]   76   17
#> [18,]   76   33
#> [19,]   76   53
#> [20,]   76    1
#> [21,]   79    3
#> [22,]   79   35
#> [23,]   79    1
#> [24,]   79    8
#> [25,]   79   11
#> [26,]   81    5
#> [27,]   81   15
#> [28,]   81   38
#> [29,]   81   11
#> [30,]   84    2
#> [31,]   84    7
#> [32,]   95    3
#> [33,]   95   51
#> [34,]  106   50
#> [35,]  124   14
#> [36,]  124   16
```
