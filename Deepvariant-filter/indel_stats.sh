#!/bin/bash
# Calculate indel statistics

if [ $# -ne 2 ]; then
    echo "Usage: $0 <indel.vcf> <out.txt>"
    exit 1
fi

indel_vcf=$1
out_txt=$2

awk 'BEGIN{n=0; total_len=0}
     !/^#/ {
         len = (length($4) > length($5)) ? length($4)-length($5) : length($5)-length($4);
         n++;
         total_len += len;
     }
     END {
         print "Indel_Count:", n;
         print "Total_Indel_Length:", total_len;
         print "Avg_Indel_Length:", (n>0 ? total_len/n : 0);
     }' "$indel_vcf" > "$out_txt"
