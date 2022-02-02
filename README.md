[![Build Status](https://travis-ci.com/damiendevienne/phylter.svg?branch=master)](https://travis-ci.com/damiendevienne/phylter)
[![R-CMD-check](https://github.com/damiendevienne/phylter/workflows/R-CMD-check/badge.svg)](https://github.com/damiendevienne/phylter/actions)
# Phylter, a tool for analyzing, visualizing and filtering phylogenomics datasets. 

**Phylter** is a tool that allows detecting, removing and visualizing outliers in phylogeneomics dataset by iteratively removing taxa in genes and optimizing a score of concordance between individual matrices.   
**Phylter** relies on Distatis (Abdi et al, 2005), an extension of multidimensional scaling to 3 dimensions to compare multiple distance matrices at once.  
**Phylter** takes as input either a collection of phylogenetic trees (that are converted to distance matrices by **Phylter**), or a collection of pairwise distance matrices (obtained from multiple sequence alignements, for instance).  
**Phylter** accepts data with missing values (missing taxa in some genes).  
**Phyler** detects outliers with a method proposed by Hubert & Vandervieren (2008) for skewed data.  
**Phylter** does not accept that the same taxa is present multiple times in the same gene. 


**Phylter** is written in R language.

## Installation
**Phylter** is not yet on CRAN. To install the development version:    


1. Install the release version of devtools from CRAN with `install.packages("devtools")`.    
2. Launch R and type:
```R
library(devtools)
install_github("damiendevienne/phylter")
```
3. Once installed the package can be loaded using:
```R
library("phylter")
```

## Usage
1. Read trees from external file and save as a list called ```trees```
```R
trees<-read.tree("treefile.tre")
```
2. (optional) Read or get gene names somewhere (same order as the trees) and save it as a vector called ```names```
3. Run phylter on your trees (see details below for possible options)
```R
results<-phylter(trees, gene.names=names)
```
>#### Options
>The phylter() function is called as follows by default: 
>```R
>phylter(X, bvalue=0, distance="patristic", k=3, k2=k, Norm=TRUE, gene.names=NULL, test.island=FALSE, verbose=TRUE, stop.criteria=1e-5)
>```
>Possible options are:    
>```bvalue```: If X is a list of trees, nodes with a support below 'bvalue' will be collapsed prior to the outlier detection.  
>```distance``` If X is a list of trees, type of distance used to compute the pairwise matrices for each tree. Can be "patristic" (sum of branch lengths separating tips, the default) or nodal (number of nodes separating tips).  
>```k``` Strength of outlier detection. The higher this value the less outliers detected.  
>```k2``` Same as k for complete gene outlier detection. To preserve complete genes from being discarded, k2 can be increased . By default, k2 = k. (see above)
>```Norm``` Should the matrices be normalized. If TRUE (the default), each matrix is divided by its median. This ensures that fast-evolving genes are not considered outliers.
>```gene.names``` List of gene names used to rename elements in X. If NULL (the default), elements are named 1,2,..,length(X).   
>```test.island``` If TRUE, only the highest value in an 'island' of outliers is considered an outlier. This prevents non-outliers hitchhiked by outliers to be considered outliers themselves (default to FALSE).   
>```verbose``` If TRUE (the default), messages are written during the filtering process to get information on what is happening  
>```stop.criteria``` The optimization stops when the gain (quality of compromise) between round *n* and round *n*+1 is smaller than this value. Default to 1e-5.  
4. Analyze the results 
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
5. Save the results
Save the results of the analysis to an external file, for example to perform cleaning on raw alignments based the results from phylter. 
```R
write.phylter(results, file="phylter.out")
```
## Example
A fungal dataset comprised of  246 genes for 21 species (Aguileta *et al.* 2008) is included in the package. To load it and test **phylter** on it: 
```R
data(fungi)
results<-phylter(fungi, distance="nodal", thres=0.2) #for example
```
   
   
---
### references
Abdi, H., Valentin, D., O’Toole, A.J., & Edelman, B. (2005). DISTATIS: The analysis of multiple distance matrices. Proceedings of the IEEE Computer Society: International Conference on Computer Vision and Pattern Recognition. (San Diego, CA, USA). pp. 42–47. https://www.utdallas.edu/~herve/abdi-distatis2005.pdf

G Aguileta, S Marthey, H Chiapello, M.-H Lebrun, F Rodolphe, E Fournier, A Gendrault-Jacquemard, T Giraud, Assessing the Performance of Single-Copy Genes for Recovering Robust Phylogenies, Systematic Biology, Volume 57, Issue 4, August 2008, Pages 613–627, https://doi.org/10.1080/10635150802306527

Hubert, M. and Vandervieren, E. (2008). An adjusted boxplot for skewed distributions. Computational Statistics and Data Analysis, 52, 5186-5201. 

de Vienne D.M., Ollier S. et Aguileta G. (2012) Phylo-MCOA: A Fast and Efficient Method to Detect Outlier Genes and Species in Phylogenomics Using Multiple Co-inertia Analysis. Molecular Biology and Evolution 29 : 1587-1598. (**This is the ancestor of Phylter**)

---
>For comments, suggestions and bug reports, please open an issue on this github repository.

