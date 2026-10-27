# Environment and installation

The tested runtime is R 4.4.3 Patched (2025-02-28 r88218) on Linux x86_64. `renv.lock` pins the package environment:

| Package | Version |
|---|---|
| AIPW | 0.6.9.3 |
| SuperLearner | 2.0-40 |
| drcmd | 0.1.0 |
| mice | 3.17.0 |
| ranger | 0.18.0 |
| rootSolve | 1.8.2.4 |
| renv | 1.3.1 |

`drcmd` is installed from `keithbarnatchez/drcmd` at commit `fd1f563867ab3827f4dba9daea8610f5c74308b7`. The lockfile also records transitive dependencies, package sources and licenses. Package versions are not recorded for every bundled replication; the lockfile specifies the tested runtime rather than certifying the environment for every result row.

## Setup

From the project root, run `Rscript scripts/setup.R`. This restores dependencies into `renv/library`, using the project-local `renv/cache` and `renv/local` directories. Installation is separate from simulation execution; libraries and caches are excluded from Git.

Setup requires access to CRAN and GitHub. Source-package installation may require a C/C++/Fortran toolchain and system libraries. HTML rendering requires Pandoc. PDF rendering additionally requires `pdflatex` and the LaTeX packages `float` and `placeins`. RStudio usually provides Pandoc; otherwise install it separately or set `RSTUDIO_PANDOC` to its directory. renv does not install these system prerequisites.

On a cluster providing the `conda_R/4.4.x` module, load it before running Rscript. Other systems can use their installed R executable.

## Scope of testing

Dependency restoration with the project library succeeds. Aggregation and HTML rendering also work from a relocated project directory with its own package library. A complete network installation on a different operating system has not been tested.

The runners check package defaults that affect model settings: AIPW outer folds, SuperLearner RF trees and drcmd internal cross-validation folds. RF-MICE passes `ntree = 500` explicitly with the ranger backend, independently of the mice package default. See [analysis settings](../report/ANALYSIS_NOTES.md). Changing package versions can change fitted values and random-number consumption even when the seed is fixed.
