pkgname <- "PhylteR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "PhylteR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('PhylteR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("PhylteR")
### * PhylteR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PhylteR
### Title: PhylteR
### Aliases: PhylteR

### ** Examples


# Detecting outliers of the dataset Fungi using nodal distances.
# This data set doesn't contain any missing data.

data(Fungi)

Results <- PhylteR(Fungi, distance = "nodal", bvalue = 0, k = 3,
thres = 0.6, gene.names = NULL, Norm = "NONE")

# See results
# Complete outliers

outgn <- Results$complete$outgn
outsp <- Results$complete$outsp

# outliers cell

outcell <- Results$CellByCell$outcell

# you can visualize the 2WR matrices (genes x species) with the function plot2WR.

plot = plot2WR(Results$Complete$mat2WR)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PhylteR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot2WR")
### * plot2WR

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot2WR
### Title: plot2WR
### Aliases: plot2WR

### ** Examples

# Detecting outliers of the dataset Fungi using nodal distances.
# This data set doesn't contain any missing data.

data(Fungi)

Results <- PhylteR(Fungi, distance = "nodal", bvalue = 0, k = 3,
thres = 0.6, gene.names = NULL, Norm = "NONE")

# you can visualize the 2WR matrices (genes x species) with the function plot2WR.

plot = plot2WR(Results$Complete$mat2WR)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot2WR", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("trees2matrices")
### * trees2matrices

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: trees2matrices
### Title: trees2matrices
### Aliases: trees2matrices

### ** Examples

# transforming a lsit of trees into a list of distances matrices using patristic distances:
data(Fungi)
matrices = trees2matrices(Fungi, distance = "patristic", bvalue = 0)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("trees2matrices", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
