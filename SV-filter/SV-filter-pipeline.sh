#!/bin/bash
# SV merging and post-analysis pipeline
# Author: Yuejingjing
# Description: Filter individual caller VCFs, merge results, adjust BND, and generate summary statistics.

# Usage:
#   bash SV-merge-pipeline.sh <depth_min> <depth_max>
#
# Example:
#   bash SV-merge-pipeline.sh 9 71

set -euo pipefail

if [ $# -ne 2 ]; then
    echo "Usage: $0 <depth_min> <depth_max>"
    exit 1
fi

depth_min=$1
depth_max=$2

echo "=== SV Merging Pipeline Started ==="
echo "Depth range: $depth_min - $depth_max"

# Step 1: Filtering
echo "[1/6] Filtering VCFs ..."

run_and_report_filter() {
    local script=$1
    local in_vcf=$2
    local out_vcf=$3
    local label=$4

    bash "$script" "$in_vcf" "$out_vcf" "$depth_min" "$depth_max" > /dev/null 2>&1

    local before=$(grep -vc "^#" "$in_vcf" || echo 0)
    local after=$(grep -vc "^#" "$out_vcf" || echo 0)
    
    local pct="0.00"
    if [ "$before" -gt 0 ]; then
        pct=$(awk -v b="$before" -v a="$after" 'BEGIN {printf "%.2f", (a/b)*100}')
    fi

    printf "  >> [%-10s] 过滤前: %5d 条 -> 过滤后: %5d 条 (保留: %s%%)\n" "$label" "$before" "$after" "$pct"
}

FILTER_BIN="/public/home/yuejingjing/tyh/protocol/SV-filter/github"

run_and_report_filter "$FILTER_BIN/filter_cutesv.sh" "cutesv.vcf" "cutesv.f.vcf" "CuteSV"
run_and_report_filter "$FILTER_BIN/filter_sniffles2.sh" "sniffles2.vcf" "sniffles2.f.vcf" "Sniffles2"
run_and_report_filter "$FILTER_BIN/filter_svim.sh" "svim.vcf" "svim.f.vcf" "SVIM"
run_and_report_filter "$FILTER_BIN/filter_pbsv.sh" "pbsv.vcf" "pbsv.f.vcf" "pbsv"

# Step 2: Merge with SURVIVOR
echo "[2/6] Merging VCFs with SURVIVOR ..."
ls SyRI.vcf svim-asm.vcf pbsv.f.vcf svim.f.vcf cutesv.f.vcf sniffles2.f.vcf > sample
SURVIVOR merge sample 1000 3 1 1 0 50 merge.vcf

# Step 3: Extract INV breakpoints
echo "[3/6] Extracting INV breakpoints ..."
bash /public/home/yuejingjing/tyh/protocol/SV-filter/github/INV-breakpoints-finder.sh merge.vcf INV.breakpoint.txt

# Step 4: Adjust BND with Python script
echo "[4/6] Adjusting BND (filter misclassified INV) ..."
python3 /public/home/yuejingjing/tyh/protocol/SV-filter/github/BND-adjust.py
grep -v -Ff deleted_ids.txt merge.vcf > merge.filtered.vcf

# Step 5: Statistics
echo "[5/6] Generating statistics ..."
bash /public/home/yuejingjing/tyh/protocol/SV-filter/github/stat_sv.sh merge.filtered.vcf statistic.number+length.txt

# Step 6: Extract DUP info
echo "[6/6] Extracting DUP info ..."
bash /public/home/yuejingjing/tyh/protocol/SV-filter/github/extract_dup.sh merge.vcf DUP.location.txt

echo "=== SV Merging Pipeline Finished Successfully ==="
echo "Results:"
echo " - merge.filtered.vcf"
echo " - statistic.number+length.txt"
echo " - DUP.location.txt"
echo " - INV.breakpoint.txt / INV.filtered.txt"
