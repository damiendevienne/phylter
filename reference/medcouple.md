# A robust measure of skewness for univariate data

Computes the medcouple, a robust measure of skewness for univariate
data. For multivariate data the medcouple is calculated on each column
of the data matrix.

## Usage

``` r
medcouple(x, do.reflect = NULL)
```

## Arguments

- x:

  An \\n\\ by \\p\\ data matrix.

- do.reflect:

  Logical indicating whether the medcouple should also be computed on
  the reflected sample `-x`, with final result
  \\(mc\\(`x`)\\-mc(-\\`x`\\))/2\\.  
  Defaults `TRUE` when \\n\<=100\\ and to `FALSE` otherwise.

## Value

mc A \\p\\-vector containing the medcouple of each column of the data
matrix `x`.

## Details

The medcouple is a robust measure of skewness yielding values between
\\-1\\ and \\1\\. For left- and right-skewed data the medcouple is
negative and positive respectively.  
The medcouple is defined as the median of the kernel function
\\h(x_i,x_j) = \frac{(x_j - med(x)) - (med(x)-x_i)}{x_j-x_i}\\ evaluated
over all couples \\(x_i,x_j)\\ where \\x_i\\ is smaller than the median
of `x` and \\x_j\\ larger than the median of `x`. When there are
multiple observations tied to the median, the kernel is defined
separately as the denominator is not defined for these observations. Let
\\m_1 \< ... \< m_k\\ denote the indices of the observations which are
tied to the median. Then \\h(x\_{m_i},x\_{m_j})\\ is defined to equal
\\-1\\ if \\i + j - 1 \< k\\, \\0\\ when \\i + j - 1 = k\\ and \\+1\\ if
\\i + j - 1 \> k\\. To compute the medcouple an algorithm with time
complexity \\O(n log(n))\\ is applied. For details, see
<https://en.wikipedia.org/wiki/Medcouple>. For numerical accuracy it is
advised, for small data sets, to compute the medcouple on both `x` and
`-x`. The final value of the medcouple may then be obtained as a linear
combination of both calculations. This procedure is warranted by the
properties of the medcouple. Indeed the medcouple of the distribution
\\X\\ equals minus the medcouple of the reflected distribution \\-X\\.
Moreover the medcouple is location and scale invariant. Note that
missing values are not allowed.

## Note

This function is extracted from the package mrfDepth - 07/2023

## References

Brys G., Hubert M., Struyf A. (2004). A robust measure of skewness.
*Journal of Computational and Graphical Statistics*, **13**, 996–1017.

## Author

P. Segaert with original code from M. Maechler and G. Brys.

## Examples

``` r
# Calculate the medcouple of univariate data sets.
# For 2000 normally distributed values
# the medcouple value is close to 0 because 
# data are not skewed
x<-rnorm(2000)
medcouple(x) 
#> [1] -0.0369128
#> attr(,"class")
#> [1] "medcouple"
# For 2000 values following a lognormal
# distribution (mean 0,sd 1), medcouple is close to 1
# because values are right-skewed
y<-rnorm(2000)
medcouple(y) 
#> [1] 0.07324025
#> attr(,"class")
#> [1] "medcouple"
# Use the option do.reflect to increase expected accuracy. 
medcouple(y, do.reflect = TRUE)
#> [1] 0.07324025
#> attr(,"class")
#> [1] "medcouple"
```
