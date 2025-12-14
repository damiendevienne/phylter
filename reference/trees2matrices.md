# Convert phylogenetic trees to distance matrices

Transform a list of phylogenetic trees into a list of phylogenetic
distance matrices.

## Usage

``` r
trees2matrices(trees, distance = "patristic", bvalue = 0)
```

## Arguments

- trees:

  A list of gene trees in format "multiphylo".

- distance:

  A method to generate distance matrices. It could be "nodal" to
  establish that the distance between two species is the number of nodes
  that separate them. Or "patristic" (default) if the distance between
  two species is the sum of branch lengths separating them. The "nodal"
  option should only be used if all species are present in all trees (no
  missing data).

- bvalue:

  This argument is only used if trees contain bootstrap values. It
  determines under what bootstrap values the nodes should be collapsed.
  Value 0 (the default) means that no nodes are collapsed.

## Value

return a list of distance matrices

## Examples

``` r
data(carnivora)
matrices<-trees2matrices(carnivora)
```
