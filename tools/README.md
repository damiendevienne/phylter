# Local audit and command-line usage

Run commands from the repository root. Temporary data and installed packages can
stay entirely in this repository:

```sh
mkdir -p .audit/tmp .audit/library .audit/downloads .audit/reference .audit/reference-lib .audit/candidate-lib
export TMPDIR="$PWD/.audit/tmp"
export R_LIBS="$PWD/.audit/library"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TZ=UTC
```

Install missing dependencies into `.audit/library`. When using `install.packages`,
use absolute `lib` and `destdir` paths; parallel installation changes directories.
Keep dependency versions identical for both installations. For reproducibility,
the audit RDS files record `sessionInfo()`.

The unchanged reference is pinned (do not substitute the working tree):

```sh
git archive 4d74241169be3882da9d5f08da47bd40604764eb | tar -x -C .audit/reference
R CMD INSTALL --library=.audit/reference-lib .audit/reference
R CMD INSTALL --library=.audit/candidate-lib .
Rscript tools/audit-run.R .audit/reference-lib .audit/reference-run 3
Rscript tools/audit-run.R .audit/candidate-lib .audit/candidate-run 3
Rscript tools/audit-compare.R .audit/reference-run.rds .audit/candidate-run.rds
```

With dependencies and the archived reference in place, `bash tools/run-audit.sh`
runs installation, numerical comparisons, kernel tests, CLI checks and separate
memory measurements. It writes logs under `.audit/` and gives each CLI test run
a fresh output directory.

`bash tools/check-package.sh` stages the candidate source inside `.audit/`, builds
it, and runs `R CMD check` with examples and both package tests. Vignettes and the
PDF manual are excluded from that check; `rmarkdown` is not needed. Staging avoids
R's attempt to recursively copy the repository into a temporary directory inside
it, while keeping all files within the repository.

Run timing processes sequentially, on the same machine, with no simultaneous
compilation. Cumulative R allocations are recorded with `Rprofmem`; they are not
peak memory and exclude native allocations. On Linux, measure peak resident
memory in separate processes:

```sh
/usr/bin/time -v Rscript tools/audit-memory.R .audit/reference-lib distatis
/usr/bin/time -v Rscript tools/audit-memory.R .audit/candidate-lib distatis
/usr/bin/time -v Rscript tools/audit-memory.R .audit/reference-lib imputation
/usr/bin/time -v Rscript tools/audit-memory.R .audit/candidate-lib imputation
```

The executable is installed with the package in its `exec` directory. R package
installation does not automatically add it to PATH. For this local installation:

```sh
export R_LIBS="$PWD/.audit/candidate-lib:$PWD/.audit/library"
export PATH="$PWD/.audit/candidate-lib/phylter/exec:$PATH"
phylter --help
phylter --trees gene_trees.nwk --out .audit/analysis
```

Alternatively run `Rscript exec/phylter ...` with that same `R_LIBS` setting.
This command calls the **installed** package; reinstall after changing R/C++
source. `--save-rds` retains the complete R result, and `--report` requests a PDF.
Neither is required for the TSV outputs.

Outputs:

| Suffix | Contents |
| --- | --- |
| `.outliers.tsv` | Headered `gene`, `species` pairs detected during optimization |
| `.discarded.tsv` | Same columns, genes discarded during initial preparation |
| `.summary.txt` | Existing human-readable package report, including both kinds of removal |
| `.session.txt` | Input path, explicit CLI options and R/dependency session information |
| `.rds` | Optional full result |
| `.pdf` | Optional plots/report |

Directory input accepts `.nwk`, `.newick`, `.tre`, `.tree`, `.treefile`, one tree
per file. Full filenames become IDs. A multi-Newick file must contain at least
two trees; parsed tree names are retained, otherwise IDs follow input order.
The parent output directory must exist. Existing outputs are refused. Patristic
distances require finite, nonnegative branch lengths. `--support` uses the input
support units directly (it does not convert proportions to percentages).

The frontend currently runs the full in-memory R analysis. It does not yet prune
trees, filter FASTA files, accept matrix files, resume interrupted runs or enforce
a memory/thread budget. Those features and the compiled-core roadmap are in
`PERFORMANCE.md`.
