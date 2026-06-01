#!/usr/bin/env bash
#SBATCH --job-name=1001_invariant_analysis
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=8:00:00
#SBATCH --account=gmonroegrp
#SBATCH --partition=bml
#SBATCH --output=logs/summarize_1001_segment_missingness_pi.%j.out
#SBATCH --error=logs/summarize_1001_segment_missingness_pi.%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/mlensink/ABCslim_trees}"
SEGMENT_SUMMARY="${SEGMENT_SUMMARY:-/group/gmonroegrp2/mlensink/ABC_data/real_data_segments/vcf_subsets_all_sites/segment_record_counts.tsv}"
OUT_DIR="${OUT_DIR:-/group/gmonroegrp2/mlensink/1001data/invariantanalysis}"
WINDOW_SIZE="${WINDOW_SIZE:-10000}"
THRESHOLDS="${THRESHOLDS:-0.50,0.75,0.80,0.90,0.95,1.00}"

cd "${PROJECT_DIR}"
mkdir -p logs "${OUT_DIR}"

python3 scripts/summarize_1001_segment_missingness_pi.py \
  --segment-summary "${SEGMENT_SUMMARY}" \
  --out-dir "${OUT_DIR}" \
  --window-size "${WINDOW_SIZE}" \
  --thresholds "${THRESHOLDS}"
