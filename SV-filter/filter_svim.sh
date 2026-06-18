#!/bin/bash
if [ $# -ne 4 ]; then
    echo "Usage: $0 <input.vcf> <output.vcf> <depth_min> <depth_max>"
    exit 1
fi

in_vcf=$1
out_vcf=$2
depth_min=$3
depth_max=$4

MIN_SUPPORT=3

awk -F'\t' -v min_dp="$depth_min" -v max_dp="$depth_max" -v min_supp="$MIN_SUPPORT" '
/^#/ { print; next }

$7 == "PASS" {
    
    supp = 0;
    if (match($8, /SUPPORT=([0-9]+)/, a)) {
        supp = a[1] + 0;
    }

    if (supp >= min_supp) {
        
        split($9, fmt, ":");
        dp_idx = 0;
        for (i=1; i<=length(fmt); i++) {
            if (fmt[i] == "DP") dp_idx = i;
        }

        if (dp_idx > 0) {
            split($10, smp, ":");
            dp = smp[dp_idx];

            if (dp == "." || (dp + 0 >= min_dp && dp + 0 <= max_dp)) {
                print $0;
            }
        }
    }
}
' "$in_vcf" > "$out_vcf"
