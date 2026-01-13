#!/bin/bash
# Split SNPs and Indels from a VCF

if [ $# -ne 3 ]; then
    echo "Usage: $0 <input.vcf> <snp.vcf> <indel.vcf>"
    exit 1
fi

in_vcf=$1
snp_vcf=$2
indel_vcf=$3

awk 'BEGIN{OFS="\t"} /^#/ {print > "'$snp_vcf'"; print > "'$indel_vcf'"; next}
     length($4)==1 && length($5)==1 {print > "'$snp_vcf'"}
     length($4)>1 || length($5)>1 {print > "'$indel_vcf'"}' "$in_vcf"

