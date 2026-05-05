#!/bin/bash
set -euo pipefail

search_dir="${1:-output}"
manifest="${2:-output/vcf_to_bgzip_manifest.txt}"

mkdir -p "$(dirname "${manifest}")"
find "${search_dir}" -type f -name '*.vcf' | sort -V > "${manifest}"

count="$(wc -l < "${manifest}")"
files_per_task="${BGZIP_FILES_PER_TASK:-25}"
echo "Wrote ${count} VCF paths to ${manifest}"
if [[ "${count}" -gt 0 ]]; then
  last_task=$(((count - 1) / files_per_task))
  echo "Submit with:"
  echo "sbatch --array=0-${last_task}%25 --export=ALL,VCF_MANIFEST=${manifest},BGZIP_FILES_PER_TASK=${files_per_task} scripts/bgzip_vcfs_array.sbatch"
fi
