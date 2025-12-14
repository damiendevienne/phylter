# print objects of class phylterfinal

Print on screen a simple description of the content of objects of class
`phylterfinal`.

## Usage

``` r
# S3 method for class 'phylterfinal'
print(x, ...)
```

## Arguments

- x:

  Object returned by function 'phylter()'.

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
print(res$Final)
#> Phylter Analysis - final state
#> List of class phylterfinal
#> 
#>    Object            Dimension
#> 1  $WR               53 x 125 
#> 2  $RV               125 x 125
#> 3  $weights          125      
#> 4  $compromise       53 x 53  
#> 5  $F                53 x 8   
#> 6  $PartialF         125      
#> 7  $species.order    53       
#> 8  $AllOptiScores    11       
#> 9  $CELLSREMOVED     94 x 2   
#> 10 $Outliers         94 x 2   
#> 11 $CompleteOutliers 2        
#> 12 $matrices         125      
#>    Content                                           
#> 1  Species x Genes reference matrix                  
#> 2  Genes x Genes RV correlation coefficients matrix  
#> 3  Weight of each gene in the compromise             
#> 4  Species x Species compromise matrix               
#> 5  Distatis coordinates of compromise                
#> 6  Distatis coordinates of gene matrices (list)      
#> 7  Name and order of species                         
#> 8  Evolution of quality of compromise (11 steps)     
#> 9  Index of cells removed (may contain imputed cells)
#> 10 Outliers detected (one row = one outlier cell)    
#> 11 Complete outliers (Gene and Species, if any)      
#> 12 Species x Species gene matrices (list)            
```
