#!/bin/bash

: "${BASE_OUT:=/home/mlensink/slimsimulations/ABCslim_trees/output}"
: "${CHROM_DIR:=chrom1}"

# Directory containing VCF files
VCF_DIR="${VCF_DIR:-${BASE_OUT}/constant/${CHROM_DIR}/vcf/}"

# Output directory for Tajima's D results
OUTPUT_DIR="${OUTPUT_DIR:-${BASE_OUT}/constant/${CHROM_DIR}/tajima}"

# Create the output directory if it doesn't already exist
mkdir -p ${OUTPUT_DIR}

# Export variables to be available in the parallel environment
export VCF_DIR OUTPUT_DIR

# Use parallel to run vcftools for each VCF file, limiting jobs to 64
parallel --jobs 64 --env VCF_DIR,OUTPUT_DIR 'vcftools --vcf {} --out ${OUTPUT_DIR}/$(basename {} .vcf) --TajimaD 100' ::: ${VCF_DIR}/*.vcf
