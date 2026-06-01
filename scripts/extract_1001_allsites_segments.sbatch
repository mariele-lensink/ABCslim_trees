#!/usr/bin/env bash
#SBATCH --job-name=1001_allsites_segments
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=6:00:00
#SBATCH --account=gmonroegrp
#SBATCH --partition=bml
#SBATCH --output=logs/extract_1001_allsites_segments.%j.out
#SBATCH --error=logs/extract_1001_allsites_segments.%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/mlensink/ABCslim_trees}"
STORAGE_DIR="${STORAGE_DIR:-/group/gmonroegrp2/mlensink/ABC_data}"

cd "${PROJECT_DIR}"

if [[ -z "${BCFTOOLS:-}" ]]; then
  if command -v bcftools >/dev/null 2>&1; then
    BCFTOOLS="$(command -v bcftools)"
  elif command -v module >/dev/null 2>&1; then
    module load bcftools/1.19
    BCFTOOLS="$(command -v bcftools)"
  elif [[ -x "/home/mlensink/anaconda3/envs/singletons/bin/bcftools" ]]; then
    BCFTOOLS="/home/mlensink/anaconda3/envs/singletons/bin/bcftools"
  else
    echo "Could not find bcftools. Set BCFTOOLS=/path/to/bcftools or load bcftools before submitting." >&2
    exit 1
  fi
fi

VCF_GZ="${VCF_GZ:-/group/gmonroegrp2/mlensink/1001data/1001genomes_snp-short-indel_with_tair10_only_ACGTN.vcf.gz}"
OUT_DIR="${OUT_DIR:-${STORAGE_DIR}/real_data_segments/vcf_subsets_all_sites}"
SUMMARY_TSV="${OUT_DIR}/segment_record_counts.tsv"

mkdir -p logs "${OUT_DIR}"

if [[ ! -s "${VCF_GZ}.csi" && ! -s "${VCF_GZ}.tbi" ]]; then
  echo "No tabix index found; building CSI index for ${VCF_GZ}"
  "${BCFTOOLS}" index --csi "${VCF_GZ}"
fi

printf "segment\tregion\trecords\tinvariant_ALT_dot\tvariant_ALT_not_dot\tvcf_gz\n" > "${SUMMARY_TSV}"

extract_region() {
  local segment="$1"
  local region="$2"
  local out_vcf="${OUT_DIR}/${segment}.vcf.gz"

  echo "Extracting ${region} -> ${out_vcf}"
  "${BCFTOOLS}" view -Oz -r "${region}" -o "${out_vcf}" "${VCF_GZ}"
  "${BCFTOOLS}" index --csi --force "${out_vcf}"

  local counts
  counts="$("${BCFTOOLS}" view -H "${out_vcf}" | awk '
    BEGIN { inv = 0; var = 0 }
    { if ($5 == ".") inv++; else var++ }
    END { print inv + var "\t" inv "\t" var }
  ')"

  printf "%s\t%s\t%s\t%s\n" "${segment}" "${region}" "${counts}" "${out_vcf}" >> "${SUMMARY_TSV}"
}

extract_region "chr1_10000000_10284120" "1:10000000-10284120"
extract_region "chr2_18000000_18280622" "2:18000000-18280622"
extract_region "chr3_22000000_22349526" "3:22000000-22349526"
extract_region "chr4_2000000_2660337" "4:2000000-2660337"
extract_region "chr5_1_377208" "5:1-377208"

echo "Wrote regional VCFs to ${OUT_DIR}"
echo "Summary: ${SUMMARY_TSV}"
