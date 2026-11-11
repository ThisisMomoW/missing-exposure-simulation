#!/usr/bin/env bash
#SBATCH --job-name=missing-exposure
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=7-00:00:00
#SBATCH --output=study-%j.out
set -euo pipefail
# Submit from the project root after making R available in the cluster environment.
args=(--run study --workers "${SLURM_CPUS_PER_TASK:-1}")
if [[ -n "${MIN_VALID_IMPUTATIONS:-}" ]]; then
  args+=(--min-valid-imputations "$MIN_VALID_IMPUTATIONS")
fi
Rscript scripts/run.R "${args[@]}"
Rscript scripts/aggregate.R --run study
Rscript scripts/figures.R --run study
