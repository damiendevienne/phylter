# Compute gene x species matrix from the result of Distatis

`Dist2WR` computes the 2WR matrix from the results obtaind with
`DistatisFast` (the fast version of distatis). For internal use mostly.

## Usage

``` r
Dist2WR(Distatis)
```

## Arguments

- Distatis:

  output of the fonction `DistatisFast`.

## Value

A matrix (gene=rows x species=col). Each cell represents a gene/species
pair, whose value represents the distance between (i) the position of
this species in this gene tree and (ii) the average position of this
species in all other gene trees.

## Examples

``` r
data(carnivora)
matrices<-phylter(carnivora, InitialOnly=TRUE, parallel=FALSE)$matrices
#> 
#> Number of Genes:    125
#> Number of Species:  53
#> --------
ds<-DistatisFast(matrices, parallel = FALSE)
WR<-Dist2WR(ds) #returns the gene x species matrix
```
