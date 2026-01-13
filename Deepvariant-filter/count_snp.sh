#!/bin/bash
# Count number of SNPs in a VCF

if [ $# -ne 2 ]; then
    echo "Usage: $0 <snp.vcf> <out.txt>"
    exit 1
fi

snp_vcf=$1
out_txt=$2

grep -v "^#" "$snp_vcf" | wc -l > "$out_txt"
