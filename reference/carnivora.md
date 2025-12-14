# 125 genes trees for 53 carnivora species

This dataset contains 125 carnivora gene trees for 53 carnivora species.
It is a subset of the dataset used in Allio et al. (2021), obtained by
randomly selecting 125 genes for which all 53 species were present.

## Usage

``` r
carnivora
```

## Format

A list of named trees of class multiPhylo

## Source

[doi:10.7554/eLife.63167](https://doi.org/10.7554/eLife.63167)

## References

Allio, R., Tilak, M. K., Scornavacca, C., Avenant, N. L., Kitchener, A.
C., Corre, E., ... & Delsuc, F. (2021). High-quality carnivoran genomes
from roadkill samples enable comparative species delineation in aardwolf
and bat-eared fox. Elife, 10, e63167.

## Examples

``` r
data(carnivora)
class(carnivora)
#> [1] "multiPhylo"
```
