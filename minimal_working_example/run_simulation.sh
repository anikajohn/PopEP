#!/usr/bin/env bash
# Simulate the minimal example end to end.
# Usage:  conda activate popep_sim && bash run_simulation.sh [threads]
set -euo pipefail
cd "$(dirname "$0")"
THREADS="${1:-4}"
CONFIG=config/scenario.yaml

python scripts/build_haplotypes.py --config "$CONFIG" --outdir results/truth
python scripts/simulate_samples.py --config "$CONFIG" --haplotypes results/truth/haplotypes.fasta \
    --outdir data --truthdir results/truth --threads "$THREADS"
python scripts/check_abundance.py --config "$CONFIG" --datadir data --truthdir results/truth \
    --outdir results/check
