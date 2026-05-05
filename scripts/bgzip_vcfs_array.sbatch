#!/bin/bash
#SBATCH -D /home/mlensink/slimsimulations/ABCslim_trees/
#SBATCH -o /home/mlensink/slimsimulations/ABCslim_trees/logs/bgzip_vcfs/%x_%A_%a.out
#SBATCH -e /home/mlensink/slimsimulations/ABCslim_trees/logs/bgzip_vcfs/%x_%A_%a.err
#SBATCH -J bgzip_vcfs
#SBATCH -t 30:00
#SBATCH --mem=750M
#SBATCH --cpus-per-task=1
#SBATCH --account=gmonroegrp
#SBATCH --partition=bmh
#SBATCH --array=0-9999%100
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=melensink@ucdavis.edu

set -euo pipefail

module load htslib/1.19.1

export VCF_MANIFEST="${VCF_MANIFEST:-/home/mlensink/slimsimulations/ABCslim_trees/output/vcf_to_bgzip_manifest.txt}"
export BGZIP_THREADS="${BGZIP_THREADS:-${SLURM_CPUS_PER_TASK}}"
export BGZIP_FILES_PER_TASK="${BGZIP_FILES_PER_TASK:-25}"

mkdir -p /home/mlensink/slimsimulations/ABCslim_trees/logs/bgzip_vcfs

if [[ ! -s "${VCF_MANIFEST}" ]]; then
  echo "Missing or empty VCF manifest: ${VCF_MANIFEST}" >&2
  exit 1
fi

total_files="$(wc -l < "${VCF_MANIFEST}")"
start_line=$((SLURM_ARRAY_TASK_ID * BGZIP_FILES_PER_TASK + 1))
end_line=$((start_line + BGZIP_FILES_PER_TASK - 1))

if (( start_line > total_files )); then
  echo "Array task ${SLURM_ARRAY_TASK_ID} starts at line ${start_line}, above ${total_files}; nothing to do."
  exit 0
fi

if (( end_line > total_files )); then
  end_line="${total_files}"
fi

echo "Compressing manifest lines ${start_line}-${end_line} of ${total_files}"
sed -n "${start_line},${end_line}p" "${VCF_MANIFEST}" | while IFS= read -r vcf_path; do
  if [[ -z "${vcf_path}" ]]; then
    continue
  fi
  if [[ ! -f "${vcf_path}" ]]; then
    echo "VCF no longer exists, skipping: ${vcf_path}"
    continue
  fi
  if [[ -f "${vcf_path}.gz" ]]; then
    echo "Compressed output already exists, skipping: ${vcf_path}.gz"
    continue
  fi

  echo "Compressing ${vcf_path} with ${BGZIP_THREADS} bgzip threads"
  bgzip -@ "${BGZIP_THREADS}" "${vcf_path}"
done
