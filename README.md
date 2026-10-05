# PhylteR, a tool for analyzing, visualizing and filtering phylogenomics datasets <img src="man/figures/logophylter.png" align="right" style="float:right; width:20%;"/>

[![CRAN_Release_Badge](https://www.r-pkg.org/badges/version-ago/phylter)](https://cran.r-project.org/package=phylter)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/phylter)](https://cran.r-project.org/package=phylter)
[![R-CMD-check](https://github.com/damiendevienne/phylter/workflows/R-CMD-check/badge.svg)](https://github.com/damiendevienne/phylter/actions)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/damiendevienne/phylter/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/damiendevienne/phylter)
[![Conda version](https://anaconda.org/damiendevienne/r-phylter/badges/version.svg)](https://anaconda.org/damiendevienne/r-phylter)



**PhylteR** detects, visualizes, and removes outlier sequences in phylogenomic
datasets. An outlier is a taxon–gene association whose evolutionary distances
are unusually discordant with the signal shared by the other genes. PhylteR can
start from a collection of gene trees or directly from pairwise distance
matrices, and it supports datasets in which some taxa are absent from some
genes.

The package is designed for a common phylogenomic problem: a sequence may be
misidentified, contaminated, paralogous, incorrectly aligned, or otherwise
inconsistent with the dominant evolutionary signal. Instead of judging each
gene tree in isolation or requiring a reference species tree, PhylteR compares
all genes jointly and identifies the particular gene–taxon cells responsible
for discordance.

## How PhylteR works

1. Gene trees are converted to patristic or nodal distance matrices. Users may
   also provide distance matrices directly through the R API.
2. Missing taxa are added and their pairwise distances are imputed from the
   other genes. Each matrix can be normalized so that differences in overall
   evolutionary rate do not dominate the comparison.
3. DISTATIS, an extension of multidimensional scaling for multiple distance
   matrices (Abdi et al. 2005), estimates a weighted compromise representing
   the signal shared across genes.
4. For every taxon in every gene, PhylteR measures the distance between its
   gene-specific position and its position in the compromise. An adjusted
   boxplot rule for skewed distributions (Hubert & Vandervieren 2008) flags
   unusually large deviations. The optional island rule avoids incorrectly
   flagging neighbouring taxa displaced by a strong outlier.
5. Candidate outliers are removed and the analysis is repeated. A proposed
   removal is retained only when it improves inter-gene concordance; iteration
   stops when the improvement falls below `stop.criteria`. Entire anomalous
   genes can also be detected using the `k2` threshold.

The result records the initial and final analyses, accepted outliers, discarded
genes, concordance scores, and the objects needed by PhylteR's summary and
visualization functions. Taxon labels must be unique within each gene.

## What is new in this version

This branch preserves PhylteR's statistical method and R API while making the
existing algorithm faster, less memory-intensive, easier to run, and easier to
validate:

- a registered C++ kernel builds the weighted DISTATIS compromise without
  allocating a complete list of weighted matrices;
- vectorized distance calculations replace repeated row-wise R callbacks;
- matrix triangle indices and normalization factors are computed once and
  reused;
- missing-distance imputation accumulates values one gene at a time instead of
  creating two additional full matrix collections;
- outlier islands are detected with a linear scan for standard taxon labels;
- two native-array leaks and two bounds-check ordering issues in the medcouple
  implementation are fixed;
- a command-line interface makes the same implementation available to users
  who do not want to write R code; and
- reference fixtures, numerical regression tests, kernel tests, CLI integration
  tests, and a reproducible performance audit now accompany the code.

On the audited Carnivora example (125 genes and 53 species), the default
analysis was **2.43× faster** (median 0.586 s versus 1.425 s). Synthetic tests
showed a **2.26×** speedup for missing-data imputation and a **39.81×** speedup
for repeated gene-to-compromise distance calculations. At 400 species × 150
genes, peak resident memory fell by 3.7% for DISTATIS and 7.2% for imputation.
These measurements were made on one documented system and should not be read as
universal performance guarantees. All reference comparisons preserved the
outlier identities, their order, the optimization trajectory, and numerical
results within tolerance. See [the performance and correctness audit](PERFORMANCE.md)
for the environment, complete results, limitations, and optimization roadmap.

The novelty of this release is therefore primarily computational and practical:
it retains the published PhylteR procedure and results while reducing avoidable
work, providing a non-interactive workflow, and adding an auditable correctness
baseline for future optimization. The current implementation still holds dense
distance matrices in memory and is not yet an out-of-core solution for very
large datasets.

PhylteR builds on Phylo-MCOA (de Vienne et al. 2012) and is implemented in R
with a small native C++ numerical core. For function documentation and a
step-by-step biological example, visit the
[PhylteR website](https://damiendevienne.github.io/phylter).

> Note: if you don't use R or don't want to use R, **containerized versions of phylter** are also available (Docker and Singularity): [https://damiendevienne.github.io/phylter/articles/phyltercontainer.html](https://damiendevienne.github.io/phylter/articles/phyltercontainer.html)

> if you use **phylter**, please cite: Comte, A., Tricou, T., Tannier, E., Joseph, J., Siberchicot, A., Penel, S., Allio, R., Delsuc, F., Dray, S., de Vienne, D.M. (2023). PhylteR: Efficient Identification of Outlier Sequences in Phylogenomic Datasets, Molecular Biology and Evolution, 40(11) msad234, [https://doi.org/10.1093/molbev/msad234](https://doi.org/10.1093/molbev/msad234)


## Installation

The current release of **PhylteR** is available on CRAN.


Installation is as easy as typing what follows at the R command prompt: 
```R
install.packages("phylter")
```

To test the optimized version described above before it reaches CRAN, install
this branch from GitHub:

1. Install the release version of `remotes` from CRAN:
```R
install.packages("remotes")
```

2. Install the development version of `phylter` from GitHub:
```R
remotes::install_github(
  "damiendevienne/phylter",
  ref = "perf/matrix-optimizations-cli-audit"
)

```
3. Once installed, the package can be loaded:
```R
library("phylter")
```

> PhylteR requires R 4.0 or later. Package installation also requires the
> system libraries needed by its R dependencies.

## Usage

A command-line frontend is included in `exec/phylter`. It accepts either a
multi-Newick file or a directory containing one Newick tree per gene:

```sh
phylter --trees gene_trees.nwk --out analysis
phylter --trees gene_trees/ --out analysis --report --save-rds
```

The CLI uses the installed R package, so its scientific results are the same as
those of the R API. It writes separate, headered TSV files for detected outliers
and genes discarded during preparation, plus a text summary and session
metadata. `--report` adds a PDF report and `--save-rds` saves the complete R
result. Existing output files are not overwritten. Run `phylter --help` for all
options, or see [local installation and CLI usage](tools/README.md) for setup and
the complete output contract.

Here is a brief introduction to the use `phylter` on a collection of gene trees. For more detailed explanations and a use case example, please visit  https://damiendevienne.github.io/phylter/.
<!-- For more more detailed examples, please go to [ADD LINK TO THE AUTOMATICALLY GENERATED WEBSITEWEB](ADD LINK TO THE AUTOMATICALLY GENERATED WEBSITEWEB). -->


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
>The `phylter` function is called as follows by default:
>```R
>phylter(X, bvalue = 0, distance = "patristic", k = 3, k2 = k, Norm = "median", 
>  Norm.cutoff = 0.001, gene.names = NULL, test.island = TRUE, 
>  verbose = TRUE, stop.criteria = 1e-5, InitialOnly = FALSE, normalizeby = "row", 
>  parallel = TRUE)
>```
>
>Arguments are as follows:
>
>- `X`: A list of phylogenetic trees (phylo object) or a list of distance matrices. Trees can have different number of leaves and matrices can have different dimensions. If this is the case, missing values are imputed.
>- `bvalue`: If `X` is a list of trees, nodes with a support below `bvalue` will be collapsed prior to the outlier detection.
>- `distance`: If `X` is a list of trees, type of distance used to compute the pairwise matrices for each tree. Can be "patristic" (sum of branch lengths separating tips, the default) or "nodal" (number of nodes separating tips).
>- `k`: Strength of outlier detection. The higher this value the less outliers detected.
>- `k2`: Same as `k` for complete gene outlier detection. To preserve complete genes from being discarded, `k2` can be increased. By default, `k2 = k`.
>- `Norm`:  Should the matrices be normalized prior to the complete analysis and how. If "median", matrices are divided by their median; if "mean", they are divided by their mean; if "none", no normalization if performed. Normalizing ensures that fast-evolving (and slow-evolving) genes are not treated as outliers. Normalization by median is a better choice as it is less sensitive to outlier values.
>- `Norm.cutoff`: Value of the median (if `Norm = "median"`) or the mean (if `Norm = "mean"`) below which matrices are simply discarded from the analysis. This prevents dividing by 0, and allows getting rid of genes that contain mostly branches of length 0 and are therefore uninformative anyway. Discarded genes, if any, are listed in the output (`out$DiscardedGenes`).
>- `gene.names`: List of gene names used to rename elements in `X`. If NULL (the default), elements are named 1,2,...,length(X).
>- `test.island`: If `TRUE` (the default), only the highest value in an *island* of outliers is considered an outlier. This prevents non-outliers hitchhiked by outliers to be considered outliers themselves.
>- `verbose`: If `TRUE` (the default), messages are written during the filtering process to get information on what is happening.
>- `stop.criteria`: The optimization stops when the gain (quality of compromise) between round *n* and round *n*+1 is smaller than this value. Default to 1e-5.
>- `InitialOnly`: Logical. If `TRUE`, only the Initial state of the data is computed.
>- `normalizeby`: Should the gene x species matrix be normalized prior to outlier detection, and how.
>- `parallel`: Logical. Should the computations be parallelized when possible? Default to `TRUE`. Note that the number of threads cannot be set by the user when `parallel = TRUE`. It uses all available cores on the machine. 

**4.** Analyze the results

To get the list of outliers detected by `phylter`, simply type:

```R
results$Final$Outliers
```

In addition, many functions allow looking at the outliers detected and comparing before and after *phy*ltering. 

```R
# Get a summary: nb of outliers, gain in concordance, etc.
summary(results)

# Show the number of species in each gene, and how many per gene are outliers
plot(results, "genes") 

# Show the number of genes where each species is found, and how many are outliers
plot(results, "species") 

# Compare before and after genes x species matrices, highlighting missing data and outliers 
# identified (not efficient for large datasets)
plot2WR(results) 

# Plot the dispersion of data before and after outlier removal. One dot represents one 
# gene x species association
plotDispersion(results) 

# Plot the genes x genes matrix showing pairwise correlation between genes
plotRV(results) 

# Plot optimization scores during optimization
plotopti(results) 
``` 

**5.** Save the results of the analysis to an external file, for example to perform cleaning on raw alignments or pruning gene trees based on the results from `phylter`.

```R
write.phylter(results, file = "phylter.out")
```


## References

- Abdi, H., O’Toole, A.J., Valentin, D. & Edelman, B. (2005). *DISTATIS: The analysis of multiple distance matrices.* Proceedings of the IEEE Computer Society: International Conference on Computer Vision and Pattern Recognition (San Diego, CA, USA). doi: 10.1109/CVPR.2005.445. https://www.utdallas.edu/~herve/abdi-distatis2005.pdf

- Allio, R., Tilak, M. K., Scornavacca, C., Avenant, N. L., Kitchener, A. C., Corre, E., ... & Delsuc, F. (2021). High-quality carnivoran genomes from roadkill samples enable comparative species delineation in aardwolf and bat-eared fox. Elife, 10, e63167. https://doi.org/10.7554/eLife.63167

- Comte, A., Tricou, T., Tannier, E., Joseph, J., Siberchicot, A., Penel, S., Allio, R., Delsuc, F., Dray, S., de Vienne, D.M. (2023). PhylteR: Efficient Identification of Outlier Sequences in Phylogenomic Datasets, Molecular Biology and Evolution, 40(11), msad234, [https://doi.org/10.1093/molbev/msad234](https://doi.org/10.1093/molbev/msad234)

- Hubert, M. and Vandervieren, E. (2008). *An adjusted boxplot for skewed distributions.* Computational Statistics and Data Analysis. https://doi.org/10.1016/j.csda.2007.11.008

- de Vienne D.M., Ollier S. et Aguileta G. (2012). *Phylo-MCOA: A Fast and Efficient Method to Detect Outlier Genes and Species in Phylogenomics Using Multiple Co-inertia Analysis.* Molecular Biology and Evolution. https://doi.org/10.1093/molbev/msr317 (This is the ancestor of phylter). 


---
For comments, suggestions and bug reports, please open an [issue](https://github.com/damiendevienne/phylter/issues) on this GitHub repository.
