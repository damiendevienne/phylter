# Phylter performance and usability audit

Reference: package 0.9.12, Git revision
`4d74241169be3882da9d5f08da47bd40604764eb`. The reference source was archived
before changes and is installed separately from the candidate. No commits or
pushes are part of this work. Local dependencies, builds, results and profiles
are under the ignored `.audit/` directory.

## Recommendation

Keep the R API and add a small C++ numerical core, with a command-line frontend
that does not require users to write R. Optimize memory layout and algorithms
before considering a complete rewrite. The current code already uses compiled
matrix products, an eigensolver and a C++ medcouple implementation.

A Python translation could improve integration with Python workflows, but it
would not by itself remove the dominant matrix storage or pairwise computation.
NumPy also delegates its numerical linear algebra to BLAS/LAPACK; that is evidence
about its architecture, not a benchmark against this package. See the
[NumPy documentation](https://numpy.org/doc/stable/reference/routines.linalg.html).
Use a shared C++ core if Python bindings eventually become necessary, rather than
maintaining two independent implementations of the scientific method.

This first implementation changes allocation patterns and a few local kernels.
It does **not** implement an out-of-core solver or claim to solve the largest
datasets. The existing eigenvalue selection, normalization definitions, island
rules, and iterative acceptance criterion are preserved.

## Measured results

Measured locally on an Intel Core i7-1165G7 with R 4.5.2, GCC/G++ 15.2,
Linux x86-64, reference BLAS
`libblas.so.3.12.1`, and `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`,
`MKL_NUM_THREADS` set to 1. Both versions use the same installed dependencies.
These are local measurements, not projected production-scale speedups.
Key dependencies were ape 5.8-1, Rfast 2.1.5.2, RSpectra 0.16-2,
Rcpp 1.1.1-1.1 and reshape2 1.4.5.

| Workload | Reference | Candidate | Speedup |
| --- | ---: | ---: | ---: |
| Full Carnivora, defaults (median of 3) | 1.425 s | 0.586 s | 2.43× |
| DISTATIS, 200 species × 80 genes (median of 3) | 0.280 s | 0.244 s | 1.15× |
| Imputation, 200 species × 80 genes, 15% missing taxa/gene (median of 3) | 0.237 s | 0.105 s | 2.26× |
| WR, 200 species × 80 genes, 20 calls (median of 3 batches) | 0.637 s | 0.016 s | 39.81× |
| Full Carnivora with missing taxa (single run) | 1.806 s | 0.873 s | 2.07× |

Default full-run samples were 1.693/1.403/1.425 s for reference and
0.777/0.586/0.582 s for candidate; first calls include warm-up. Small kernel
timings, especially the 16 ms WR batch, have limited precision. Other configuration
timings in `.audit/comparison.log` are single observations, not robust estimates.

In separate processes at **400 species × 150 genes**, maximum resident memory
(including R, dependencies and input generation) changed as follows:

| Workload | Reference peak RSS | Candidate peak RSS | Reduction |
| --- | ---: | ---: | ---: |
| DISTATIS | 979,516 KiB | 942,788 KiB | 3.7% |
| Imputation | 933,380 KiB | 866,584 KiB | 7.2% |

At 200 × 80, cumulative R allocation for DISTATIS fell from 304,770,792 to
197,333,048 bytes (35.3%). Imputation allocation **increased** from 249,425,904
to 263,021,712 bytes (5.5%), even though its smaller live working set reduced
peak RSS in the separate larger test. Cumulative allocations and peak memory
measure different things; these changes are not a solution to the dense-storage
limit. `Rprofmem` also does not track all native allocations.

Sampling the original full Carnivora run attributed about 53% of time to outlier
detection and 32% to DISTATIS (inclusive, overlapping call costs). After the linear
island scan, those shares were about 19% and 73%. The matrix algebra now dominates
this example. Larger `N`/`K` can have very different bottlenecks.

**Correctness:** all 12 complete/initial-only configurations and synthetic kernels
passed old/new comparisons. Default output remains **94 outliers, 11 accepted
states**, with concordance approximately **0.8623535 → 0.9442600**. Comparisons
check identities/order and numerical states, not only these summary numbers.
The CLI integration checks passed for file/directory inputs, names containing
spaces, TSV outputs, PDF generation, invalid arguments and overwrite protection.
`tests/fixtures/carnivora-reference.rds` keeps the original result for future CI.
The final candidate also passes **R CMD check: Status OK**, including examples,
`kernels.R` and `regression.R`. The PDF manual and vignette rebuild were excluded;
the suggested `rmarkdown` dependency is unavailable in this local environment.

## Complete source inventory

All 15 R files, all four original native/build files, package metadata and
namespace, README, both vignettes, generated manual pages, and CI configuration
were inspected. Existing example output/report files are documentation artifacts,
not an executable regression suite. There were no automated tests in the source
at the reference revision.

| Source | Role and finding |
| --- | --- |
| `R/phylter.R` | Orchestrates preparation, DISTATIS, clustering, detection and repeated candidate evaluation. Retains initial and final large objects. Replaces all previously accepted outliers on each proposal, not just new ones. Computes full candidate projections even when quality will reject them. Clusters even when island detection is disabled and during the whole-gene pass. |
| `R/DistatisFast.R` | Main dense algebra. Centers every matrix, packs triangles, builds a gene-by-gene Gram/RV matrix, finds weights, builds a compromise, selects axes, projects every gene. `auto` asks for `N-1` eigenpairs, despite documentation emphasizing partial eigenanalysis. |
| `R/Dist2WR.R` | Computes per-species distances between gene and compromise projections. The original invokes an R function separately for every row of every gene. |
| `R/PreparePhylterData.R` | Converts trees, imputes/reorders, normalizes and discards low-scale genes. Computes every median/mean twice, normalizes before discarding, and retains original matrices. |
| `R/impMean.R` | Expands every gene to the union of taxa. Original summation materializes full lists of zero-filled matrices and presence masks. Has an unreliable fallback for taxa never observed together; see correctness findings. |
| `R/detect.outliers.R` | Uses a global adjusted boxplot after optional row/column normalization. Island detection constructs all pairs of flagged positions just to find adjacent runs, then merges lists. Repeated `rbind` adds allocations. Whole-gene detection uses low weights. |
| `R/normalize.R` | Uses `apply` and medians. Zero medians can produce nonfinite values. Dimension dropping matters for small inputs. Naming of species/genes options is easy to misinterpret. |
| `R/trees2matrices.R` | Uses `ape` for support collapse and cophenetic matrices. Builds a second tree list before conversion. Support collapse acts on edges whose parent is the low-support node; biological intent should be checked. |
| `R/medcouple.R` | Validates/converts input and calls native code through `.C`. Default reflection contradicts its documentation. |
| `src/medcouple.cpp` | Thin `.C` bridge; includes Eigen although this file uses no Eigen calculations. |
| `src/mc.cpp` | Actual medcouple selection algorithm. Already avoids materializing all pairs. Original early returns leak two arrays. Many unused Eigen-related declarations/includes remain from imported code. |
| `src/phylter_init.c`, `src/Makevars` | Native registration and BLAS/LAPACK linking. Original registration declares an unused external function. Several compilation dependencies merit a separate minimal-build audit. |
| `R/summary.phylter.R` | Summary and four print methods. No-outlier outputs need robust zero/empty handling; static print descriptions can disagree with actual field order. Reported percentage gain is a percentage-point difference, not relative percent improvement. |
| `R/plot.phylter.R` | Six plotting functions. Expands matrices into data frames and performs potentially expensive clustering. `plot2WR` clusters even with `clust=FALSE`, and uses a gene reordering for missing-data markers even without reordering the heatmap. Large heatmaps need aggregation/rasterization. |
| `R/write.phylter.R` | Human-readable text and optional PDF. Useful foundation but not a clean machine-output contract. PDF device cleanup is not protected on errors; empty PDF filename fallback is assigned but not used. |
| `R/rename.genes.R` | Assigns names without checking count, uniqueness or empty labels. Preparation only calls it when names are absent, so explicit replacement names are ignored for already named inputs. |
| `R/simtrees.R` | Development simulator. Repeated list growth and an unbounded topology-changing loop; small and degenerate inputs need validation. Some simulated whole-gene/species information is not carried into returned outlier truth. |
| `R/data.R`, `data/carnivora.rda` | Included biological example: 125 genes, 53 species. Useful reference, but too small and too complete to establish large-data or imputation correctness alone. |
| Documentation / metadata / CI | Container vignette already describes a Python wrapper and pruning/filtering scripts, but those scripts are absent here. Integrate with that interface deliberately. CI checks the R package; add regression/CLI coverage. README's blanket claim that R requires GSL should be replaced by specific dependency requirements. |

## Where time and memory go

Let `K` be genes, `N` the union of species, `r` retained axes, and `T` the number
of accepted plus rejected candidate evaluations.

| Work | Approximate time | Important storage |
| --- | --- | --- |
| Expand, impute, normalize, center | `O(K N²)` per pass (median implementation also matters) | Each dense collection costs `8 K N²` bytes |
| Full pairwise gene similarities | `O(K² N²)` per evaluation | RV alone costs `8 K²` bytes; packed features about `4 K N²` |
| Leading RV eigenpair | Iterative, roughly `O(iterations × K²)` after constructing RV | Solver workspace plus RV |
| Compromise accumulation | `O(K N²)` | Original weighted collection adds `8 K N²`; now one output matrix |
| Automatic compromise spectrum | Near-full spectrum, cubic-scale work in `N` | `O(N²)` vectors/workspace |
| Gene projection and WR | `O(K N² r)` | `O(K N r)` partial projections |
| Detection / islands | Medcouple plus quantiles; original island pair construction quadratic in flagged positions per gene | Pair tables, masks and output lists |

Repeated evaluation multiplies much of this work by `T`. A collection with
1,000 genes and 1,000 species is **7.45 GiB for one dense double collection**.
10,000 genes and 500 species require **18.63 GiB for one collection**, plus
**0.75 GiB for one RV matrix**. These are storage arithmetic, not measured peak
RAM estimates. Multiple live collections, triangle buffers, eigensolver workspace,
projections and R allocation/garbage collection increase the peak. R objects can
share storage until modified, so simply counting variable names overestimates
some copies and misses temporary allocations.

Tree distance and imputed matrices are generally dense: converting them to a
sparse format is unlikely to help. More threads can increase both workspace and
memory bandwidth contention. Benchmark thread counts rather than using every
core by default.

## Implemented first pass

* `Dist2WR`: `sqrt(rowSums(...))` replaces row-wise `apply` callbacks.
* Triangle indices are computed once per DISTATIS evaluation, not once per gene.
* A registered C++ `.Call` kernel accumulates the weighted compromise directly
  into one newly allocated matrix, retaining gene order and leaving inputs intact.
* Preparation reuses normalization factors already computed.
* Named, distinct taxon labels use a linear scan of consecutive outlier positions
  instead of pair tables and nested island merging. Legacy behavior is retained
  for unusual labels (including the old `"out"` sentinel).
* Imputation sums values and presence counts one gene at a time, avoiding two
  temporary full collections. The expanded input/output collections still exist.
* Medcouple's two early-allocated work arrays now have automatic lifetime; bounds
  are checked before reading indexed values in two loops.
* `verbose=FALSE` no longer emits the unconditional final stop message.
* An executable `exec/phylter` frontend accepts Newick files/directories,
  validates common input errors and writes separate tables, a summary and session
  metadata. It calls the same package implementation, not a reimplementation.

The standalone allocation-balance test in `tools/medcouple-sanitizer.cpp` reports
1,320 unreleased raw arrays for the reference and zero for the candidate on the
same repeated small/constant/general inputs. AddressSanitizer and UBSan pass for
the candidate with leak detection disabled. LeakSanitizer itself cannot run under
this environment's process tracing; the allocation-balance test supplies the
direct old/new lifetime comparison instead.

## Next optimizations, in priority order

1. **Split candidate evaluation into quality and projection stages.** Compute
   centered matrices, RV and its leading eigenpair first. Only compute the
   compromise spectrum and partial projections after a proposal is accepted.
   Preserve the existing gain test and output. This saves rejected evaluations
   without changing the method, subject to numerical trajectory checks.
2. **Native streaming/fused matrix kernels.** Center/pack in one pass; preallocate
   feature buffers instead of `lapply` plus `cbind`; stream imputation directly
   from original smaller matrices. Avoid retaining full centered collections when
   only projections/WR are required. Keep the median over the full matrix,
   including diagonal zeros and duplicated symmetric entries, for compatibility.
3. **Finish detection cleanup.** The linear island scan is now implemented and
   preserves group ordering; 8,190 exhaustive small patterns and 200 randomized/
   sentinel cases match the original helper. Still remove repeated result `rbind`
   and skip clustering when island detection is disabled or only gene weights
   are used. Preserve all tied maxima and exact `CELLSREMOVED` ordering.
4. **Optional compact results.** A CLI filtering run needs outlier IDs, discarded
   IDs, scores, summaries and provenance. It need not retain both dense states,
   RV matrices and every projection. Add an explicit retention option while
   keeping the existing rich R return value as the compatibility default. Plotting
   should clearly report which retained data it needs.
5. **Matrix-free RV for many genes.** For centered symmetric matrix `S_g`, form a
   vector of its diagonal and `sqrt(2)` times its upper triangle, normalized by
   the Frobenius norm. Stack those vectors as columns of `Z`. Then `RV = ZᵀZ` and
   `RV v = Zᵀ(Z v)`. The leading eigenpair can be computed without storing `K²`
   entries or explicitly computing every pairwise similarity. Each operator
   application costs `O(K N²)`, with convergence-dependent repetition.
   [RSpectra already accepts matrix-vector functions](https://spectralib.org/r-interface).
   A native operator avoids repeated R allocation. Full RV output must be optional.
6. **Packed/block storage and a memory budget.** Packed symmetric matrices roughly
   halve persistent matrix storage. Blocked feature products bound temporary
   memory; memory mapping/out-of-core iteration trades RAM for I/O. Account for
   original missing-taxon masks and replacements. Merely changing the file format
   while still expanding all matrices in RAM will not solve the problem.
7. **Cache unaffected genes.** Centering, norms and unchanged RV blocks can be
   reused, but the affected set includes every gene with a previously removed
   cell: its replacement changes with the compromise each round. Caching based
   only on newly detected cells is incorrect. RV weights and the compromise still
   change globally.
8. **Revisit automatic axes separately.** A fixed small rank or randomized solver
   can save substantial work but can change the outlier list. The current
   broken-stick rule normalizes using the computed `N-1` eigenvalues and requests
   largest magnitude eigenpairs. Negative eigenvalues and near-zero modes matter;
   do not silently replace this with a positive-spectrum or fixed-rank rule.
   Test dense LAPACK versus near-full RSpectra as a compatible optimization first.

GPU support and distributed processing should come after a contiguous blocked
core and representative measurements. Repeated transfers, modest matrix sizes,
global medcouple/detection and peak device memory can erase their benefit.

## Language and distribution options

| Option | Benefit | Cost / recommendation |
| --- | --- | --- |
| Optimized R + small C++ kernels | Reuses `ape`, existing outputs and plotting; reduces copies and callback overhead | Best near-term route; already demonstrated by the compromise kernel |
| Python + NumPy/SciPy | Python pipeline integration, familiar packaging for some users | Requires matching normalization, quantiles, eigenspaces, medcouple, imputation, tree handling and iteration; no automatic complexity reduction |
| Standalone C++ executable + shared core | Precise memory ownership, predictable thread budget, R/Python bindings can share one method | Best possible long-term core if scale warrants it; substantial I/O, numerical and packaging work |
| Rscript CLI + packaged environment | Immediately usable without writing R | R runtime remains a dependency; package/executable installation must be documented clearly |

For further native bindings, [Rcpp provides R/C++ integration](https://www.rcpp.org/pdf/Rcpp-package.pdf).
The small accumulator added here uses R's existing C API and needs no new package
dependency. Do not remove Eigen/RcppArmadillo/Rfast dependencies until package
installation and all supported paths are checked on Linux, macOS and Windows.

## Biologist and pipeline user experience

The new frontend supports, for example:

```sh
phylter --trees genes.nwk --out analysis
phylter --trees gene_trees/ --out analysis --k 3 --report
```

See `tools/README.md` for local installation and verification. Gene files are
sorted deterministically; full filenames are IDs, avoiding ambiguous underscore
or dot splitting. A multi-Newick file uses its parsed tree names when available,
otherwise input-order IDs. Output files are never overwritten. The default
disables the package's uncontrolled Rfast parallel branch; `--parallel` is
explicit. This is **not** a guarantee of one BLAS thread.

For a polished tool release, add:

* `phylter check`, `phylter run`, and a separate `phylter filter` command. A check
  should report taxa/gene counts, missingness, duplicate IDs, invalid lengths,
  support-label scale, degenerate genes and estimated memory before computation.
* A real `--threads N` budget coordinated across BLAS and native loops, plus
  `--memory`, a reproducible seed/solver configuration and progress by phase.
* An explicit manifest connecting gene IDs to tree/alignment paths. Pruned trees
  and filtered FASTA files should be new outputs; preserve sequence headers and
  input alignments. Specify minimum remaining taxa and whole-gene removal rules.
* Headered TSVs with stable columns and a versioned JSON run manifest containing
  input hashes, package/core versions, command parameters and stopping reason.
* Optional reports that aggregate large datasets, with full numeric tables kept
  separately. Report generation must not be required for filtering to succeed.
* Installation through a maintained Conda/Bioconda recipe and pinned containers,
  with example Snakemake/Nextflow invocations. These are distribution proposals,
  not claims that new packages or integrations have been published.

## Correctness findings that must not be hidden inside optimization

These are source-level findings, not a claim that every path has been exhaustively
reproduced. Existing behavior remains unchanged except for the memory lifetime and
quiet-output fixes listed above.

* `ImputeFromClosestNeighbors` calls `mean(a, b)`: in R, the second positional
  argument is `trim`, not a second observation. The intended average presumably
  needs `mean(c(a, b))`, but that changes scientific output and needs a separate
  decision/fixture. Nearest-neighbor searches also need a termination rule for
  disconnected taxon co-occurrence graphs.
* `impMean` decides whether to impute from matrix dimensions. Full-size matrices
  containing explicit NA entries bypass imputation. Validate missingness explicitly.
* Medcouple documentation says reflection defaults on for `n <= 100`; code
  actually turns it on for `n > 100`. Decide whether documentation or code expresses
  intended behavior before changing it.
* Fewer than two surviving genes, very few taxa, zero matrix norms, zero retained
  eigenvalues and zero row/column medians are not handled consistently. The
  full pipeline needs an explicit supported-domain contract.
* Public matrix entry points need checks for square/symmetric/finite matrices,
  aligned unique row/column names and nonnegative distances. Unknown options
  should fail early. New CLI validation covers common tree input problems but
  does not replace validation in the R API.
* No-outlier results can be NULL; summary counts and plotting should consistently
  handle empty tables. Tiny-data shape dropping can turn a matrix into a vector.
* Scaling and imputation are scientifically ordered: imputation currently happens
  before per-gene normalization and even before removing low-scale genes. Moving
  it later would change values and must be treated as a method change.

## Verification contract

`tools/audit-run.R` runs separately installed reference and candidate packages in
separate processes using the same R and dependency versions. `audit-compare.R`
requires identical outlier IDs/order, removed-cell order, complete/discarded gene
results and species order. It compares every accepted optimization score, initial
and final WR/RV matrices, weights, compromise and stored matrices at tolerance
`1e-8`. Projection geometry is compared using Gram products, because eigenvectors
can change sign or rotate within a tied eigenspace without changing the analysis.

Coverage includes the full supplied dataset with defaults, nodal distances, mean
normalization, no normalization, no islands, a different whole-gene threshold,
missing taxa, an added discarded zero-distance gene, and initial-only mode.
The parallel matrix-product path and column/disabled WR normalization are also
compared.
Additional 200-species/80-gene synthetic kernels test imputation, DISTATIS and WR.
`tests/kernels.R` checks independent identities and native early-return paths.

Future method changes must carry two tests: comparison with reference behavior
where compatibility is promised, and explicit expected results for corrected
behavior where the old version is wrong. Matching a historical bug is not proof
of scientific correctness. Large-scale tests should separately vary taxa, genes,
missingness, outlier density and thread count. A small biological example alone
does not validate scalability, pathological eigenvalues or disconnected taxa.
