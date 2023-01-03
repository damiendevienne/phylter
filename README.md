[![Build Status](https://travis-ci.com/damiendevienne/phylter.svg?branch=master)](https://travis-ci.com/damiendevienne/phylter)
[![R-CMD-check](https://github.com/damiendevienne/phylter/workflows/R-CMD-check/badge.svg)](https://github.com/damiendevienne/phylter/actions)
# phylter, a tool for analyzing, visualizing and filtering phylogenomics datasets

**phylter** is a tool that allows detecting, removing and visualizing outliers in phylogenomics dataset by iteratively removing taxa in genes and optimizing a score of concordance between individual matrices.   
**phylter** relies on Distatis (Abdi et al, 2005), an extension of multidimensional scaling to 3 dimensions to compare multiple distance matrices at once.   
**phylter** takes as input either a collection of phylogenetic trees (that are converted to distance matrices by `phylter`), or a collection of pairwise distance matrices (obtained from multiple sequence alignements, for instance).  
**phylter** accepts data with missing values (missing taxa in some genes).  
**phylter** detects outliers with a method proposed by Hubert & Vandervieren (2008) for skewed data.  
**phylter** does not accept that the same taxa is present multiple times in the same gene. 

**phylter** is written in R language.


## Installation
**phylter** is not yet on CRAN (deposit in process). To install the development version:    

1. Install the release version of `remotes` from CRAN:
```R
install.packages("remotes")
```

2. Install the development version of `phylter` from GitHub:
```R
remotes::install_github("damiendevienne/phylter")

```
3. Once installed, the package can be loaded:
```R
library("phylter")
```

> Note: Phylter requires R version > 4.0, otherwise it cannot be installed. Also, R uses the GNU Scientific Library. On Ubuntu, this can be installed prior to the installation of the Phylter package by typing `sudo apt install libgsl-dev` in a terminal. 


## Usage
**1.** With the `read.tree` function from the `ape` package, read trees from external file and save as a list called `trees`.
```R
if (!requireNamespace("ape", quietly = TRUE))
   install.packages("ape")
trees <- ape::read.tree("treefile.tre")
```

**2.** (optional) Read or get gene names somewhere (same order as the trees) and save it as a vector called `names`.

**3.** Run `phylter` on your trees (see details below for possible options).
```R
results <- phylter(trees, gene.names = names)

```
>#### Options
>The phylter() function is called as follows by default:
>```R
>phylter(X, bvalue = 0, distance = "patristic", k = 3, k2 = k, Norm = "median", 
>  Norm.cutoff = 0.001, gene.names = NULL, test.island = TRUE, 
>  verbose = TRUE, stop.criteria = 1e-5, InitialOnly = FALSE, normalizeby = "row", parallel = TRUE)
>```
>
>Arguments are as follows:
>- `X`: A list of phylogenetic trees (phylo object) or a list of distance matrices. Trees can have different number of leaves and matrices can have different dimensions. If this is the case, missing values are imputed.
>- `bvalue`: If X is a list of trees, nodes with a support below `bvalue` will be collapsed prior to the outlier detection.
>- `distance`: If X is a list of trees, type of distance used to compute the pairwise matrices for each tree. Can be "patristic" (sum of branch lengths separating tips, the default) or "nodal" (number of nodes separating tips).
>- `k`: Strength of outlier detection. The higher this value the less outliers detected.
>- `k2`: Same as `k` for complete gene outlier detection. To preserve complete genes from being discarded, `k2` can be increased. By default, `k2 = k`.
>- `Norm`:  Should the matrices be normalized prior to the complete analysis and how. If "median", matrices are divided by their median; if "mean" they are divided by their mean; if "none", no normalization if performed. Normalizing ensures that fast-evolving (and slow-evolving) genes are not treated as outliers. Normalization by median is a better choice as it is less sensitive to outlier values.
>- `Norm.cutoff`: Value of the median (if `Norm = "median"`) or the mean (if `Norm = "mean"`) below which matrices are simply discarded from the analysis. This prevents dividing by 0, and allows getting rid of genes that contain mostly branches of length 0 and are therefore uninformative anyway. Discarded genes, if any, are listed in the output (`out$DiscardedGenes`).
>- `gene.names`: List of gene names used to rename elements in `X`. If NULL (the default), elements are named 1,2,...,length(X).
>- `test.island`: If TRUE (the default), only the highest value in an *island* of outliers is considered an outlier. This prevents non-outliers hitchhiked by outliers to be considered outliers themselves.
>- `verbose`: If TRUE (the default), messages are written during the filtering process to get information on what is happening.
>- `stop.criteria`: The optimization stops when the gain (quality of compromise) between round *n* and round *n*+1 is smaller than this value. Default to 1e-5.
>- `InitialOnly`: Logical. If TRUE, only the Initial state of the data is computed.
>- `normalizeby`: Should the gene x species matrix be normalized prior to outlier detection, and how.
>- `parallel`: Logical. Should the computations be parallelized when possible? Default to TRUE. Note that the number of threads cannot be set by the user when `parallel=TRUE`. It uses all available cores on the machine. 

**4.** Analyze the results

Many functions allow looking at the outliers detected and comparing before and after:

```R
summary(results) # Get a summary: nb of outliers, gain in concordance, etc.
plot(results, "genes") # show the number of species in each gene, and how many per gene are outliers
plot(results, "species") # show the number of genes where each species is found, and how many are outliers
plot2WR(results) # compare before and after genes x species matrices, highlighting missing data and outliers identified. (not efficient for large datasets)
plotDispersion(results) # plot the dispersion of data before and after outlier removal. One dot represents one gene x species association.
plotRV(results) # plot the genes x genes matrix showing pairwise correlation between genes.
plotopti(results) #plot optimization scores during optimization.
``` 

**5.** Save the results

Save the results of the analysis to an external file, for example to perform cleaning on raw alignments based on the results from `phylter`.
```R
write.phylter(results, file = "phylter.out")
```


## Example
A `carnivora` dataset comprised of 125 genes for 53 species (Allio et al. 2021) is included in the package. To load it and test `phylter` on it: 
```R
data(carnivora)
results <- phylter(carnivora) # for example
```
   
   
---
## References

- Abdi, H., O’Toole, A.J., Valentin, D. & Edelman, B. (2005). *DISTATIS: The analysis of multiple distance matrices.* Proceedings of the IEEE Computer Society: International Conference on Computer Vision and Pattern Recognition (San Diego, CA, USA). doi: 10.1109/CVPR.2005.445. https://www.utdallas.edu/~herve/abdi-distatis2005.pdf

- Hubert, M. and Vandervieren, E. (2008). *An adjusted boxplot for skewed distributions.* Computational Statistics and Data Analysis. doi: 10.1016/j.csda.2007.11.008. 

- de Vienne D.M., Ollier S. et Aguileta G. (2012). *Phylo-MCOA: A Fast and Efficient Method to Detect Outlier Genes and Species in Phylogenomics Using Multiple Co-inertia Analysis.* Molecular Biology and Evolution. doi: 10.1093/molbev/msr317 (This is the ancestor of phylter).


---
For comments, suggestions and bug reports, please open an [issue](https://github.com/damiendevienne/phylter/issues) on this GitHub repository.

