#!/usr/bin/env bash
set -euo pipefail
# Dependencies must already be installed; invoke from the repository root.
export TMPDIR="$PWD/.audit/tmp"
export R_LIBS="$PWD/.audit/library"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export TZ=UTC
mkdir -p .audit/reference-lib .audit/candidate-lib
R CMD INSTALL --library=.audit/reference-lib .audit/reference > .audit/reference-install.log 2>&1
R CMD INSTALL --library=.audit/candidate-lib . > .audit/candidate-install.log 2>&1
Rscript tools/audit-run.R .audit/reference-lib .audit/reference-run 3 > .audit/reference-run.log 2>&1
Rscript tools/audit-run.R .audit/candidate-lib .audit/candidate-run 3 > .audit/candidate-run.log 2>&1
Rscript tools/audit-compare.R .audit/reference-run.rds .audit/candidate-run.rds > .audit/comparison.log 2>&1
R_LIBS="$PWD/.audit/candidate-lib:$R_LIBS" Rscript tests/kernels.R > .audit/kernels.log 2>&1
R_LIBS="$PWD/.audit/candidate-lib:$R_LIBS" Rscript tests/regression.R > .audit/regression.log 2>&1
Rscript tools/test-islands.R .audit/reference > .audit/islands.log 2>&1
cli_dir=$(mktemp -d "$PWD/.audit/cli.XXXXXX")
Rscript tools/test-cli.R .audit/candidate-lib .audit/reference-run.rds "$cli_dir" > .audit/cli.log 2>&1
for mode in distatis imputation; do
  for version in reference candidate; do
    /usr/bin/time -v -o ".audit/$version-$mode-memory.log" \
      Rscript tools/audit-memory.R ".audit/$version-lib" "$mode" \
      > ".audit/$version-$mode-memory.stdout" 2>&1
  done
done
cat .audit/comparison.log .audit/cli.log
