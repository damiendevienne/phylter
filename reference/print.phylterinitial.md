# print objects of class phylterinitial

Print on screen a simple description of the content of objects of class
`phylterinitial`.

## Usage

``` r
# S3 method for class 'phylterinitial'
print(x, ...)
```

## Arguments

- x:

  Object present in the `$initial` element of the object returned by
  function [`phylter()`](phylter.md).

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
print(res$Initial)
#> Phylter Analysis - initial state
#> List of class phylterinitial
#> 
#>   Object      Dimension Content                                         
#> 1 $mat.data   125       List of original distance matrices, one per gene
#> 2 $WR         53 x 125  Species x Genes reference matrix                
#> 3 $RV         125 x 125 Genes x Genes RV correlation coefficients matrix
#> 4 $weights    125       Weight of each gene in the compromise           
#> 5 $compromise 53 x 53   Species x Species compromise matrix             
#> 6 $F          53 x 6    Distatis coordinates of compromise              
#> 7 $matrices   125       Distatis coordinates of gene matrices (list)    
#> 8 $PartialF   125       Species x Species gene matrices (list)          
```
