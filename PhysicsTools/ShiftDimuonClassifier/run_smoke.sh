#!/usr/bin/env bash
# Standalone read-only SM pilot; no scram/cmsRun/Condor or shared-runtime edits.
set -eo pipefail
source /cvmfs/sft.cern.ch/lcg/views/LCG_108/x86_64-el9-gcc13-opt/setup.sh
set -u
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
study_directory=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd "$study_directory"
dataset_directory=${1:-artifacts/sm_v10_40files_local}
results_directory=${2:-artifacts/pilot_seed42}
python3 -m unittest discover -s tests -v
python3 train.py --dataset "$dataset_directory" --output "$results_directory" --epochs 100 --seed 42
python3 validate_data_path.py --dataset "$dataset_directory" --models "$results_directory/models" \
  --output "$results_directory/truth_removed"
python3 plot_results.py --dataset "$dataset_directory" --results "$results_directory"
