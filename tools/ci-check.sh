#!/bin/bash
set -euo pipefail
mkdir -p ci-results
rpm -qa | sort > ci-results/rpm-versions.txt
icpx --version > ci-results/compiler.txt
R --version > ci-results/r-version.txt
Rscript tools/check-docs.R
chmod +x configure
R CMD build . > ci-results/build.log 2>&1
R CMD INSTALL sbyadanear_0.4.0.tar.gz > ci-results/install.log 2>&1
lib=$(Rscript -e 'cat(system.file("libs", "sbyadanear.so", package="sbyadanear"))')
ldd "$lib" > ci-results/linkage.txt
grep -q libiomp5 ci-results/linkage.txt
grep -q libmkl_intel_thread ci-results/linkage.txt
if grep -q libgomp ci-results/linkage.txt; then
  echo 'Unexpected GNU OpenMP linkage' >&2; exit 1
fi
Rscript -e 'library(sbyadanear); print(sby_hpc_cpu_report()); testthat::test_dir("tests/testthat", reporter="summary", stop_on_failure=TRUE)' \
  > ci-results/tests.log 2>&1
Rscript tools/thread-evidence.R > ci-results/threads.log 2>&1
python3 tools/check-mkl-verbose.py ci-results/threads.log
R CMD check --no-manual sbyadanear_0.4.0.tar.gz > ci-results/check.log 2>&1
cp sbyadanear.Rcheck/00check.log ci-results/00check.log
if grep -Eq '[1-9][0-9]* (WARNING|ERROR)' ci-results/00check.log; then exit 1; fi
