#!/usr/bin/env bash
set -euo pipefail
# Stage source beside the temporary directory: R CMD build cannot copy a source
# tree into its own descendant when TMPDIR is inside that tree.
root="$PWD"
export TMPDIR="$root/.audit/tmp" R_LIBS="$root/.audit/library"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 TZ=UTC
export _R_CHECK_FORCE_SUGGESTS_=false
stage=$(mktemp -d "$root/.audit/check-source.XXXXXX")
mkdir -p .audit/candidate-build
tar --exclude='./.audit' --exclude='./.git' --exclude='./.aider*' \
    --exclude='./.agents' --exclude='./.codex' --exclude='./src/*.o' \
    --exclude='./src/*.so' -cf - . | tar -xf - -C "$stage"
cd .audit/candidate-build
R CMD build --no-build-vignettes --no-manual "$stage" > build.log 2>&1
R CMD check --no-manual --no-vignettes --ignore-vignettes phylter_0.9.12.tar.gz > check.log 2>&1
cat phylter.Rcheck/00check.log
