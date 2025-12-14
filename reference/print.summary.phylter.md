# print summary of phylter objects

Prints on screen the summary of an object of class `summmary.phylter` as
returned by the [`summary.phylter()`](summary.phylter.md) function.

## Usage

``` r
# S3 method for class 'summary.phylter'
print(x, ...)
```

## Arguments

- x:

  Object returned by function 'summary.phylter()'.

- ...:

  Additional arguments.

## Value

NA

## Examples

``` r
data(carnivora)
res <- phylter(carnivora, parallel = FALSE)
#> 
#> Number of Genes:    125
#> Number of Species:  53
#> --------
#> Initial score: 0.86235
#>     28  new cells to remove -> New score: 0.90272 -> OK
#>     18  new cells to remove -> New score: 0.90833 -> OK
#>     16  new cells to remove -> New score: 0.91501 -> OK
#>     18  new cells to remove -> New score: 0.92561 -> OK
#>     5  new cells to remove -> New score: 0.93404 -> OK
#>     4  new cells to remove -> New score: 0.93692 -> OK
#>     2  new cells to remove -> New score: 0.93712 -> OK
#>     1  new cells to remove -> New score: 0.94392 -> OK
#>     1  new cells to remove -> New score: 0.94417 -> OK
#>     1  new cells to remove -> New score: 0.94426 -> OK
#>  => No more outliers detected  ->  Checking for complete gene outliers
#>  => No more outliers detected  ->  STOPPING OPTIMIZATION
#> --------
#> 
#> Total number of outliers detected: 94
#>   Number of complete gene outliers : 0
#>   Number of complete species outliers : 0
#> 
#> Gain (concordance between matrices): 8.19% 
#> Loss (data filtering): 1.42% 
summary <- summary(res)
print(summary)
#> 
#> Total number of outliers detected: 94
#>   Number of complete gene outliers : 0
#>   Number of complete species outliers : 0
#> 
#> Gain (concordance between matrices): 8.19% 
#> Loss (data filtering): 1.42% 
```
