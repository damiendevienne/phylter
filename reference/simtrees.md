# Simplistic simulation of gene trees with outliers

Simple (simplistic) simulation of trees with outliers.

## Usage

``` r
simtrees(Ngn, Nsp, Nsp.out = 0, Ngn.out = 0, Nb.cell.outlier = 0, 
brlen.sd = 0, out.type="topology",bl.mult=2)
```

## Arguments

- Ngn:

  Number of gene trees to simulate.

- Nsp:

  Number of species (tips) per tree.

- Nsp.out:

  Number of outlier species (also called rogue taxa). 0 = none.

- Ngn.out:

  Number of outlier genes. 0 = none.

- Nb.cell.outlier:

  Number of times one species in one tree is an outlier. 0 = none. The
  type of outlier is set by `out.type=`.

- brlen.sd:

  Heterogeneity of branch lengths in trees. A value with mean 0 and
  standard deviation equal to brlen.sd is added to each branch length.
  If the resulting branch length has a negative value, its absolute
  value is taken.

- out.type:

  The type of cell outlier. Can be "topology" (the default), where
  outlier species are simulated by branching them elsewhere (randomly)
  in the tree, "brlength", where terminal branches of outlier species
  are multiplied by `bl.mult` (see after) or "both" where half of the
  outliers are "topology" and the other half are "brlength". In the
  latter case, if the number of outliers (`Nb.cell.outlier`) id odd,
  there is one more topology outlier than brlength outlier simulated.

- bl.mult:

  Multiplier of terminal branches of outlier species when
  `out.type="topology"` or `out.type="both"`. Ignored otherwise.

## Value

A list X containing a list of trees in `multiPhylo` format (X\$trees)
and a list of outliers (X\$outl).

## Details

The simulation process is as follows: a first tree is generated with the
`rtree()` function and is then duplicated and modified according to the
parameters chosen by the user.

## Examples

``` r
 
# Very basic simulator, for debugging purpose mainly.
# examples: 30 genes, 120 species, 2 outlier species, 3 outlier genes
# 4 gene/species outliers of type "topology", branch length variance = 0.6.
# The branch length multiplier is set to 2 but this is
# ignored if out.type="topology"
simu<-simtrees(30,120,2,3,4,0.6, "topology",2)
trees<-simu$trees #list of trees
outl<-simu$outl #list of simulated outliers and their type


```
