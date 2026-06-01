#!/usr/bin/env bash
#SBATCH --job-name=1001_allsites_het
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=8:00:00
#SBATCH --account=gmonroegrp
#SBATCH --partition=bmh
#SBATCH --array=0-4
#SBATCH --output=logs/heterozygosity_1001_allsites/%x_%A_%a.out
#SBATCH --error=logs/heterozygosity_1001_allsites/%x_%A_%a.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-/home/mlensink/ABCslim_trees}"
VCF_DIR="${VCF_DIR:-/group/gmonroegrp2/mlensink/ABC_data/real_data_segments/vcf_subsets_all_sites}"
OUT_DIR="${OUT_DIR:-/group/gmonroegrp2/mlensink/ABC_data/real_data_segments/heterozygosity}"
TMP_ROOT="${TMP_ROOT:-${SLURM_TMPDIR:-/tmp}}"

cd "${PROJECT_DIR}"
mkdir -p logs/heterozygosity_1001_allsites "${OUT_DIR}"

if [[ -z "${VCFTOOLS_BIN:-}" ]]; then
  if command -v vcftools >/dev/null 2>&1; then
    VCFTOOLS_BIN="$(command -v vcftools)"
  elif command -v module >/dev/null 2>&1; then
    module load vcftools/0.1.16
    VCFTOOLS_BIN="$(command -v vcftools)"
  else
    echo "Could not find vcftools. Set VCFTOOLS_BIN=/path/to/vcftools before submitting." >&2
    exit 1
  fi
fi

if [[ -z "${BCFTOOLS_BIN:-}" ]]; then
  if command -v bcftools >/dev/null 2>&1; then
    BCFTOOLS_BIN="$(command -v bcftools)"
  elif command -v module >/dev/null 2>&1; then
    module load bcftools/1.19 || true
    BCFTOOLS_BIN="$(command -v bcftools || true)"
  else
    BCFTOOLS_BIN=""
  fi
fi

segments=(
  "chr1_10000000_10284120"
  "chr2_18000000_18280622"
  "chr3_22000000_22349526"
  "chr4_2000000_2660337"
  "chr5_1_377208"
)

task_id="${SLURM_ARRAY_TASK_ID:-0}"
if (( task_id < 0 || task_id >= ${#segments[@]} )); then
  echo "Invalid SLURM_ARRAY_TASK_ID=${task_id}; expected 0-4" >&2
  exit 1
fi

SEGMENT="${segments[$task_id]}"
VCF="${VCF_DIR}/${SEGMENT}.vcf.gz"
TMP_DIR="${TMP_ROOT}/vcftools_het_${SLURM_JOB_ID:-manual}_${task_id}"
VCFTOOLS_OUT="${TMP_DIR}/${SEGMENT}"
SEGMENT_OUT="${OUT_DIR}/${SEGMENT}.heterozygosity.tsv"
COMBINED_OUT="${OUT_DIR}/all_segments.heterozygosity.tsv"

if [[ ! -s "${VCF}" ]]; then
  echo "Missing VCF: ${VCF}" >&2
  exit 1
fi

rm -rf "${TMP_DIR}"
mkdir -p "${TMP_DIR}"

echo "Segment: ${SEGMENT}"
echo "VCF: ${VCF}"
echo "Temporary files: ${TMP_DIR}"
echo "Output: ${SEGMENT_OUT}"
echo "vcftools: ${VCFTOOLS_BIN}"

"${VCFTOOLS_BIN}" --gzvcf "${VCF}" --012 --out "${VCFTOOLS_OUT}"

if [[ -n "${BCFTOOLS_BIN}" ]]; then
  read -r total_vcf_sites variant_alt_not_dot < <(
    "${BCFTOOLS_BIN}" view -H "${VCF}" |
      awk 'BEGIN { total = 0; variant = 0 } { total++; if ($5 != ".") variant++ } END { print total, variant }'
  )
else
  read -r total_vcf_sites variant_alt_not_dot < <(
    gzip -cd "${VCF}" |
      awk 'BEGIN { total = 0; variant = 0 } /^#/ { next } { total++; if ($5 != ".") variant++ } END { print total, variant }'
  )
fi

awk -v segment="${SEGMENT}" \
    -v total_vcf_sites="${total_vcf_sites}" \
    -v variant_alt_not_dot="${variant_alt_not_dot}" '
  {
    for (i = 2; i <= NF; i++) {
      gt = $i
      if (gt != -1) {
        called[i]++
      }
      if (gt == 1) {
        het[i]++
        alt_seen[i] = 1
      } else if (gt == 2) {
        alt_seen[i] = 1
      }
    }
  }
  END {
    snp_sites = 0
    het_sites = 0
    called_snp_genotypes = 0
    het_snp_genotypes = 0
    called_sites = 0

    for (i in called) {
      if (called[i] > 0) {
        called_sites++
      }
    }

    for (i in alt_seen) {
      snp_sites++
      called_snp_genotypes += called[i]
      het_snp_genotypes += het[i]
      if (het[i] > 0) {
        het_sites++
      }
    }

    prop_het_genotypes_at_snp_sites = (called_snp_genotypes > 0) ? het_snp_genotypes / called_snp_genotypes : "NA"
    prop_het_sites_per_total_sites = (total_vcf_sites > 0) ? het_sites / total_vcf_sites : "NA"
    prop_het_sites_per_snp_sites = (snp_sites > 0) ? het_sites / snp_sites : "NA"

    OFS = "\t"
    print "segment",
      "total_vcf_sites",
      "variant_alt_not_dot",
      "called_sites_in_012",
      "snp_sites_with_alt_genotype",
      "het_sites",
      "called_snp_genotypes",
      "het_snp_genotypes",
      "prop_het_genotypes_at_snp_sites",
      "prop_het_sites_per_total_sites",
      "prop_het_sites_per_snp_sites"
    print segment,
      total_vcf_sites,
      variant_alt_not_dot,
      called_sites,
      snp_sites,
      het_sites,
      called_snp_genotypes,
      het_snp_genotypes,
      prop_het_genotypes_at_snp_sites,
      prop_het_sites_per_total_sites,
      prop_het_sites_per_snp_sites
  }
' "${VCFTOOLS_OUT}.012" > "${SEGMENT_OUT}"

{
  flock 9
  if [[ ! -s "${COMBINED_OUT}" ]]; then
    head -n 1 "${SEGMENT_OUT}" > "${COMBINED_OUT}"
  fi
  tail -n 1 "${SEGMENT_OUT}" >> "${COMBINED_OUT}"
} 9>"${COMBINED_OUT}.lock"

rm -rf "${TMP_DIR}"

echo "Wrote ${SEGMENT_OUT}"
echo "Updated ${COMBINED_OUT}"
