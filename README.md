# Phylter, a tool for analyzing, visualizing and filtering phylogenomics datasets. 

**Phylter** is a tool that allows detecting, removing and visualizing outliers in phylogeneomics dataset by iteratively removing taxa in genes and optimizing a score of concordance between individual matrices. 
**Phylter** relies on Distatis (Abdi et al, 2005), an extension of multidimensional scaling to 3 dimensions to compare multiple distance matrices at once.
**Phylter** takes as input either a collection of phyloegenetic trees (that are converted to distance matrices by **Phylter**), or a directly a collection of pairwise distance matrices (that one can obtain from multiple sequence alignements).
**Phylter** accepts data with missing values (missing taxa in some genes). 
**Phylter** does not accept that the same taxa is present multiple times in the same gene. 


**Phylter** is written in R language.

## Installation
**Phylter** is not yet on CRAN. To use it: 
In the console: 
```bash
git clone https://github.com/damiendevienne/phylter.git
cd phylter
R
```
And within R:
```R
# Load packages (install beforehand if not installed)
require(ape)
require(RSpectra)
require(ggplot2)
require(reshape2)

# source functions
source("R/trees2matrices.R")
source("R/rename.genes.R")
source("R/impMean.R")
source("R/DistatisFast.R")
source("R/Dist2WR.R")
source("R/normalize.R")
source("R/detect.outliers.R")
source("R/phylter.R")
source("R/simulate.trees.R")
source("R/summary.phylter.R")
source("R/plot2WR.R")

# Ready to use phylter!
```

## Usage
1. Read trees from external file and save as a list called ```trees```
2. (optional) Read or get gene names somewhere (same order as the trees) and save it as a vector called ```names```
3. Run phylter on your trees (see details below for possible options) by typing```results<-phylter(trees, gene.names=names)```
4. Analyze the results with ```summary(results)``` and all plotting functions available (```plot(results)```, ```plot2WR(results)```, ```plotDispersion(results)```, ```plotRV(results)```, ```plotopti(results)```, )
5. Write the results to an external file for further analysis and raw data cleaning with ```writeOutput(results)```

## Options of the ```phylter()``` function
The phylter() function is called as follows by default: 
```r
phylter(X, bvalue=0, distance="patristic", k=3, thres=0.3, Norm=TRUE, keep.species=TRUE, gene.names=NULL, test.island=TRUE, verbose=TRUE, stop.criteria=1e-5)
```
Possible options are: 
```bvalue```: If X is a list of trees, nodes with a support below 'bvalue' will be collapsed prior to the outlier detection.
```distance``` If X is a list of trees, type of distance used to compute the pairwise matrices for each tree. Can be "patristic" (sum of branch lengths separating tips, the default) or nodal (number of nodes separating tips).
```k``` Strength of outlier detection. The higher this value the less outliers detected.
```thres``` For the detection of complete outliers. Threshold above which genes or species are considered as complete outliers. 0.3 means that a gene (resp. species) is a complete outlier if it is detected as outlier for more than 30% of its species (resp. genes).
```Norm``` Should the matrices be normalized. If TRUE (the default), each matrix is normalized such that its first eigenvalie is equal to one.
```keep.species``` If TRUE, species are protected from being detected as complete outliers and filtered out. 
```gene.names``` List of gene names used to rename elements in X. If NULL (the default), elements are named 1,2,..,length(X). 
```test.island``` This should not be modified. If TRUE (the default), only the highest value in an 'island' of outliers is considered an outlier. This prevents non-outliers hitchhiked by outliers to be considered outliers themselves. 
```verbose``` If TRUE (the default), messages are written during the filtering process to get information on what is happening
```stop.criteria``` The optimization stops when the gain (quality of compromise) between round *n* and round *n*+1 is smaller than this value. Default to 1e-5.






---
### references
Abdi, H., Valentin, D., O’Toole, A.J., & Edelman, B. (2005). DISTATIS: The analysis of multiple distance matrices. Proceedings of the IEEE Computer Society: International Conference on Computer Vision and Pattern Recognition. (San Diego, CA, USA). pp. 42–47.


