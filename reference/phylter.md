# Filter phylogenomics datasets

Detection and filtering out of outliers in a list of trees or a list of
distance matrices.

## Usage

``` r
phylter(
  X,
  bvalue = 0,
  distance = "patristic",
  k = 3,
  k2 = k,
  Norm = "median",
  Norm.cutoff = 0.001,
  gene.names = NULL,
  test.island = TRUE,
  verbose = TRUE,
  stop.criteria = 1e-05,
  InitialOnly = FALSE,
  normalizeby = "row",
  parallel = TRUE
)
```

## Arguments

- X:

  A list of phylogenetic trees (phylo object) or a list of distance
  matrices. Trees can have different number of leaves and matrices can
  have different dimensions. If this is the case, missing values are
  imputed.

- bvalue:

  If X is a list of trees, nodes with a support below 'bvalue' will be
  collapsed prior to the outlier detection.

- distance:

  If X is a list of trees, type of distance used to compute the pairwise
  matrices for each tree. Can be "patristic" (sum of branch lengths
  separating tips, the default) or nodal (number of nodes separating
  tips). The "nodal" option should only be used if all species are
  present in all genes.

- k:

  Strength of outlier detection. The higher this value the less outliers
  detected (see details).

- k2:

  Same as k for complete gene outlier detection. To preserve complete
  genes from being discarded, k2 can be increased. By default, k2 = k
  (see above).

- Norm:

  Should the matrices be normalized prior to the complete analysis and
  how. If "median" (the default), matrices are divided by their median,
  if "mean" they are divided by their mean, if "none", no normalization
  if performed. Normalizing ensures that fast-evolving (and
  slow-evolving) genes are not treated as outliers. Normalization by
  median is a better choice as it is less sensitive to outlier values.

- Norm.cutoff:

  Value of the median (if `Norm="median"`) or the mean (if
  `Norm="mean"`) of phylogenetic distance matrices below which genes are
  simply discarded from the analysis. This prevents dividing by 0, and
  allows getting rid of genes that contain mostly branches of length 0
  and are therefore uninformative anyway. Discarded genes, if any, are
  listed in the output `out$DiscardedGenes`.

- gene.names:

  List of gene names used to rename elements in X. If NULL (the
  default), 0 elements are named 1,2,..,length(X).

- test.island:

  If TRUE (the default), only the highest value in an 'island' of
  outliers is considered an outlier. This prevents non-outliers
  hitchhiked by outliers to be considered outliers themselves.

- verbose:

  If TRUE (the default), messages are written during the filtering
  process to get information of what is happening

- stop.criteria:

  The optimisation stops when the gain in concordance between matrices
  between round `n` and round `n+1` is smaller than this value. Default
  to 1e-5.

- InitialOnly:

  Logical. If TRUE, only the Initial state of the data is computed.

- normalizeby:

  Should the gene x species matrix be normalized prior to outlier
  detection, and how.

- parallel:

  Should the computations be parallelized when possible? Default to
  TRUE. Note that the number of threads cannot be set by the user when
  \`parallel=TRUE\`. It uses all available cores on the machine.

## Value

A list of class 'phylter' with the 'Initial' (before filtering) and
'Final' (after filtering) states, or a list of class 'phylterinitial'
only, if InitialOnly=TRUE. The function also returns the list of
DiscardedGenes, if any.

## Examples

``` r
data(carnivora)

# using default paramaters
res <- phylter(carnivora, parallel = FALSE) # perform the phylter analysis
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
res # brief summary of the analysis
#> Phylter Analysis
#> List of class phylter
#> 
#> Call: phylter(X = carnivora, parallel = FALSE)
#> 
#> $Initial Initial matrices and values, before optimization
#> $Final       Final matrices, scores, outliers, after optimization
#> $DiscardedGenes  List of discarded genes (not analyzed by phylter)
#> 
#> 
#> 
#> Tips:
#>    Use summary(x) to get an overview of the results.
#>    Use plot(x) to see the distribution of outliers
#>    Use plot2WR(x) to compare WR matrices before and after
#>    Use write.phylter(x) to write the results to an easily parsable file
#>    Use plotDispersion(x) to compare Distatis projections before and after
res$DiscardedGenes # list of genes discarded prior to the analysis
#> character(0)
res$Initial # See all elements prior to the analysis
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
res$Final # See all elements at the end of the analysis
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
res$Final$Outliers # Print all outliers detected
#>       [,1]                      [,2]                        
#>  [1,] "ENSG00000005381_MPO"     "Arctocephalus_gazella"     
#>  [2,] "ENSG00000005381_MPO"     "Ursus_maritimus"           
#>  [3,] "ENSG00000005381_MPO"     "Ailurus_fulgens"           
#>  [4,] "ENSG00000106511_MEOX2"   "Panthera_tigris"           
#>  [5,] "ENSG00000106511_MEOX2"   "Mustela_putorius"          
#>  [6,] "ENSG00000114686_MRPL3"   "Procyon_lotor"             
#>  [7,] "ENSG00000116157_GPX7"    "Vulpes_vulpes"             
#>  [8,] "ENSG00000116761_CTH"     "Otocyon_megalotis"         
#>  [9,] "ENSG00000120053_GOT1"    "Phoca_vitulina"            
#> [10,] "ENSG00000120053_GOT1"    "Gulo_gulo"                 
#> [11,] "ENSG00000123307_NEUROD4" "Arctocephalus_gazella"     
#> [12,] "ENSG00000132254_ARFIP2"  "Panthera_leo"              
#> [13,] "ENSG00000132254_ARFIP2"  "Odobenus_rosmarus"         
#> [14,] "ENSG00000132254_ARFIP2"  "Felis_catus"               
#> [15,] "ENSG00000132254_ARFIP2"  "Arctocephalus_gazella"     
#> [16,] "ENSG00000132693_CRP"     "Hyaena_hyaena"             
#> [17,] "ENSG00000132693_CRP"     "Pteronura_brasiliensis"    
#> [18,] "ENSG00000132693_CRP"     "Neovison_vison"            
#> [19,] "ENSG00000132693_CRP"     "Eumetopias_jubatus"        
#> [20,] "ENSG00000133135_RNF128"  "Cryptoprocta_ferox"        
#> [21,] "ENSG00000133135_RNF128"  "Panthera_tigris"           
#> [22,] "ENSG00000133135_RNF128"  "Phoca_vitulina"            
#> [23,] "ENSG00000134240_HMGCS2"  "Gulo_gulo"                 
#> [24,] "ENSG00000134240_HMGCS2"  "Spilogale_gracilis"        
#> [25,] "ENSG00000138675_FGF5"    "Taxidea_taxus"             
#> [26,] "ENSG00000138675_FGF5"    "Lutra_lutra"               
#> [27,] "ENSG00000143125_PROK1"   "Potos_flavus"              
#> [28,] "ENSG00000149573_MPZL2"   "Paradoxurus_hermaphroditus"
#> [29,] "ENSG00000073111_MCM2"    "Paradoxurus_hermaphroditus"
#> [30,] "ENSG00000106511_MEOX2"   "Manis_javanica"            
#> [31,] "ENSG00000106511_MEOX2"   "Acinonyx_jubatus"          
#> [32,] "ENSG00000106511_MEOX2"   "Panthera_pardus"           
#> [33,] "ENSG00000106511_MEOX2"   "Leptonychotes_weddellii"   
#> [34,] "ENSG00000114686_MRPL3"   "Arctocephalus_gazella"     
#> [35,] "ENSG00000116157_GPX7"    "Canis_familiaris"          
#> [36,] "ENSG00000116157_GPX7"    "Otocyon_megalotis"         
#> [37,] "ENSG00000116761_CTH"     "Vulpes_vulpes"             
#> [38,] "ENSG00000132693_CRP"     "Enhydra_lutris"            
#> [39,] "ENSG00000132693_CRP"     "Crocuta_Crocuta"           
#> [40,] "ENSG00000132693_CRP"     "Arctocephalus_gazella"     
#> [41,] "ENSG00000132693_CRP"     "Zalophus_californianus"    
#> [42,] "ENSG00000133135_RNF128"  "Ursus_americanus"          
#> [43,] "ENSG00000138675_FGF5"    "Mustela_putorius"          
#> [44,] "ENSG00000138675_FGF5"    "Mellivora_capensis"        
#> [45,] "ENSG00000143125_PROK1"   "Paradoxurus_hermaphroditus"
#> [46,] "ENSG00000149573_MPZL2"   "Mungos_mungo"              
#> [47,] "ENSG00000106511_MEOX2"   "Panthera_onca"             
#> [48,] "ENSG00000106511_MEOX2"   "Felis_catus"               
#> [49,] "ENSG00000106511_MEOX2"   "Enhydra_lutris"            
#> [50,] "ENSG00000106511_MEOX2"   "Ailuropoda_melanoleuca"    
#> [51,] "ENSG00000106511_MEOX2"   "Odobenus_rosmarus"         
#> [52,] "ENSG00000106511_MEOX2"   "Paradoxurus_hermaphroditus"
#> [53,] "ENSG00000106511_MEOX2"   "Neofelis_nebulosa"         
#> [54,] "ENSG00000114686_MRPL3"   "Taxidea_taxus"             
#> [55,] "ENSG00000116157_GPX7"    "Lycaon_pictus"             
#> [56,] "ENSG00000132693_CRP"     "Lutra_lutra"               
#> [57,] "ENSG00000133135_RNF128"  "Procyon_lotor"             
#> [58,] "ENSG00000138675_FGF5"    "Enhydra_lutris"            
#> [59,] "ENSG00000143125_PROK1"   "Nasua_narica"              
#> [60,] "ENSG00000144355_DLX1"    "Acinonyx_jubatus"          
#> [61,] "ENSG00000149573_MPZL2"   "Suricata_suricatta"        
#> [62,] "ENSG00000149573_MPZL2"   "Helogale_parvula"          
#> [63,] "ENSG00000106511_MEOX2"   "Lynx_pardinus"             
#> [64,] "ENSG00000106511_MEOX2"   "Gulo_gulo"                 
#> [65,] "ENSG00000106511_MEOX2"   "Ursus_maritimus"           
#> [66,] "ENSG00000106511_MEOX2"   "Eumetopias_jubatus"        
#> [67,] "ENSG00000106511_MEOX2"   "Neomonachus_schauinslandi" 
#> [68,] "ENSG00000106511_MEOX2"   "Potos_flavus"              
#> [69,] "ENSG00000106511_MEOX2"   "Panthera_leo"              
#> [70,] "ENSG00000106511_MEOX2"   "Prionailurus_bengalensis"  
#> [71,] "ENSG00000106511_MEOX2"   "Canis_familiaris"          
#> [72,] "ENSG00000106511_MEOX2"   "Neovison_vison"            
#> [73,] "ENSG00000114686_MRPL3"   "Pteronura_brasiliensis"    
#> [74,] "ENSG00000116761_CTH"     "Canis_familiaris"          
#> [75,] "ENSG00000133135_RNF128"  "Mellivora_capensis"        
#> [76,] "ENSG00000133135_RNF128"  "Acinonyx_jubatus"          
#> [77,] "ENSG00000138675_FGF5"    "Gulo_gulo"                 
#> [78,] "ENSG00000138675_FGF5"    "Neovison_vison"            
#> [79,] "ENSG00000143125_PROK1"   "Gulo_gulo"                 
#> [80,] "ENSG00000149573_MPZL2"   "Cryptoprocta_ferox"        
#> [81,] "ENSG00000106511_MEOX2"   "Lynx_canadensis"           
#> [82,] "ENSG00000106511_MEOX2"   "Zalophus_californianus"    
#> [83,] "ENSG00000106511_MEOX2"   "Mirounga_angustirostris"   
#> [84,] "ENSG00000133135_RNF128"  "Pteronura_brasiliensis"    
#> [85,] "ENSG00000143125_PROK1"   "Otocyon_megalotis"         
#> [86,] "ENSG00000106511_MEOX2"   "Puma_concolor"             
#> [87,] "ENSG00000106511_MEOX2"   "Callorhinus_ursinus"       
#> [88,] "ENSG00000143125_PROK1"   "Vulpes_vulpes"             
#> [89,] "ENSG00000143196_DPT"     "Ailurus_fulgens"           
#> [90,] "ENSG00000138675_FGF5"    "Pteronura_brasiliensis"    
#> [91,] "ENSG00000143125_PROK1"   "Lycaon_pictus"             
#> [92,] "ENSG00000143125_PROK1"   "Manis_javanica"            
#> [93,] "ENSG00000116761_CTH"     "Lycaon_pictus"             
#> [94,] "ENSG00000143196_DPT"     "Taxidea_taxus"             

# \donttest{
# Change the call to phylter to use nodal distances, instead of patristic: 
res <- phylter(carnivora, distance = "nodal")
#> 
#> Number of Genes:    125
#> Number of Species:  53
#> --------
#> Initial score: 0.90793
#>     4  new cells to remove -> New score: 0.90917 -> OK
#>  => No more outliers detected  ->  Checking for complete gene outliers
#>  => No more outliers detected  ->  STOPPING OPTIMIZATION
#> --------
#> 
#> Total number of outliers detected: 4
#>   Number of complete gene outliers : 0
#>   Number of complete species outliers : 0
#> 
#> Gain (concordance between matrices): 0.12% 
#> Loss (data filtering): 0.06% 
# }
```
