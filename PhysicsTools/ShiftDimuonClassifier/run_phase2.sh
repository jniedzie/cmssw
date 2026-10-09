#!/usr/bin/env bash
# Reproduce the fixed second comparison on the completed fresh SM table.
set -eo pipefail
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh
set -u
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/jniedzie/shift_classifier_mpl
study_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$study_directory"
dataset_directory=${1:-artifacts/sm_v10_phase2_fresh}
results_directory=${2:-artifacts/phase2_seed71_72}
python3 -m unittest discover -s tests -v
python3 phase2_train.py --dataset "$dataset_directory" --output "$results_directory" \
  --epochs 120 --seeds 71 72 --split-seed 71
python3 validate_data_path.py --dataset "$dataset_directory" --models "$results_directory/models" \
  --output "$results_directory/truth_removed"
python3 acceptance_audit.py --dataset "$dataset_directory" --results "$results_directory" \
  --output "$results_directory/acceptance_audit"
python3 plot_results.py --dataset "$dataset_directory" --results "$results_directory"
