#!/usr/bin/env bash
set -euo pipefail

VCF_GZ="${VCF_GZ:-/home/mlensink/rawdata/arabidopsis/1001genomes_snp-short-indel_only_ACGTN.vcf.gz}"
GFF_DIR="${GFF_DIR:-/home/mlensink/slimsimulations/ABCslim_trees/genomeinfo}"
OUT_BASE="${OUT_BASE:-/home/mlensink/slimsimulations/ABCslim_trees/output/real_data_segments}"
TAJIMA_WINDOW="${TAJIMA_WINDOW:-100}"

SUBSET_DIR="${OUT_BASE}/vcf_subsets"
TAJIMA_DIR="${OUT_BASE}/tajima"
SUMMARY_TSV="${OUT_BASE}/segment_snp_counts.tsv"

mkdir -p "${SUBSET_DIR}" "${TAJIMA_DIR}"

gff_files=(
  "${GFF_DIR}/TAIR10_GFF3_genes_subset_gene100_chr1_10000000.gff"
  "${GFF_DIR}/TAIR10_GFF3_genes_subset_gene100_chr2_18000000.gff"
  "${GFF_DIR}/TAIR10_GFF3_genes_subset_gene100_chr3_22000000.gff"
  "${GFF_DIR}/TAIR10_GFF3_genes_subset_gene100_chr4_2000000.gff"
  "${GFF_DIR}/TAIR10_GFF3_genes_subset_gene100.gff"
)

parse_region() {
  local gff_file="$1"
  awk '
    /^##subset / {
      chrom = start = end = ""
      for (i = 1; i <= NF; i++) {
        split($i, pair, "=")
        if (pair[1] == "chrom") chrom = pair[2]
        if (pair[1] == "start") start = pair[2]
        if (pair[1] == "end") end = pair[2]
      }
      if (chrom != "" && start != "" && end != "") {
        print chrom "\t" start "\t" end
        exit
      }
    }
  ' "${gff_file}"
}

normalize_chrom_for_vcf() {
  local chrom="$1"
  chrom="${chrom#Chr}"
  chrom="${chrom#chr}"
  printf "%s\n" "${chrom}"
}

printf "segment\tchrom\tstart\tend\ttotal_snp_count\ttajima_file\n" > "${SUMMARY_TSV}"

for gff_file in "${gff_files[@]}"; do
  if [[ ! -f "${gff_file}" ]]; then
    echo "Missing GFF: ${gff_file}" >&2
    exit 1
  fi

  region_info="$(parse_region "${gff_file}")"
  if [[ -z "${region_info}" ]]; then
    echo "Could not parse subset header from ${gff_file}" >&2
    exit 1
  fi

  IFS=$'\t' read -r chrom start end <<< "${region_info}"
  vcf_chrom="$(normalize_chrom_for_vcf "${chrom}")"
  segment="$(basename "${gff_file}" .gff)"
  subset_prefix="${SUBSET_DIR}/${segment}"
  tajima_prefix="${TAJIMA_DIR}/${segment}"
  subset_vcf="${subset_prefix}.recode.vcf"
  tajima_file="${tajima_prefix}.Tajima.D"

  echo "Subsetting ${chrom}:${start}-${end} as ${vcf_chrom}:${start}-${end} -> ${subset_vcf}"
  tabix -h "${VCF_GZ}" "${vcf_chrom}:${start}-${end}" > "${subset_vcf}"

  snp_count="$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "${subset_vcf}")"
  if [[ "${snp_count}" -gt 0 ]]; then
    echo "Calculating Tajima's D for ${segment}"
    vcftools \
      --vcf "${subset_vcf}" \
      --out "${tajima_prefix}" \
      --TajimaD "${TAJIMA_WINDOW}"
  else
    printf "CHROM\tBIN_START\tN_SNPS\tTajimaD\n" > "${tajima_file}"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${segment}" "${chrom}" "${start}" "${end}" "${snp_count}" "${tajima_file}" >> "${SUMMARY_TSV}"
done

echo "Wrote outputs to ${OUT_BASE}"
echo "Summary: ${SUMMARY_TSV}"
