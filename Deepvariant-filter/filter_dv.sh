#!/bin/bash
# Filter VCF by DP field

if [ $# -ne 4 ]; then
    echo "Usage: $0 <input.vcf.gz> <output.vcf> <minDP> <maxDP>"
    exit 1
fi

in_vcf=$1
out_vcf=$2
depth_min=$3
depth_max=$4

zcat "$in_vcf" | awk -F'\t' -v min=$depth_min -v max=$depth_max 'BEGIN{OFS="\t"}
/^#/ {print; next}
{
    split($9,f,":"); split($10,s,":");
    for(i in f){
        if(f[i]=="DP" && s[i]!="." && s[i]+0 >= min && s[i]+0 <= max){
            print; next
        }
    }
}' > "$out_vcf"
