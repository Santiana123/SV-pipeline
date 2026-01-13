#!/bin/bash
# SNP/Indel filtering and statistics pipeline
# Author: Yuejingjing
# Description: Filter DeepVariant VCF, split SNP/Indel, and generate statistics.

# Usage:
#   bash SNP-Indel-pipeline.sh <sample_name> <depth_min> <depth_max>
#
# Example:
#   bash SNP-Indel-pipeline.sh Kapoho 7 53

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <sample_name> <depth_min> <depth_max>"
    exit 1
fi

name=$1
depth_min=$2
depth_max=$3

echo "=== SNP/Indel Pipeline Started ==="
echo "Sample: $name"
echo "Depth range: $depth_min - $depth_max"

# 输入文件
in_vcf="${name}.deepvariant.g.vcf.gz"

# 输出文件
filtered_vcf="${name}.filtered.vcf"
snp_vcf="${name}.snp.vcf"
indel_vcf="${name}.indel.vcf"
snp_count="${name}.snp.count.txt"
indel_stats="${name}.indel.statistic.txt"


# Step 1: DP过滤
echo "[1/4] Filtering by depth ..."
bash /public/home/stu_tongyihan/project/7.solo_and_wild_vs_sunset/scripts/filter_dv.sh \
    "$in_vcf" "$filtered_vcf" $depth_min $depth_max

# Step 2: 拆分SNP/Indel
echo "[2/4] Splitting SNPs and Indels ..."
bash /public/home/stu_tongyihan/project/7.solo_and_wild_vs_sunset/scripts/split_snp_indel.sh \
    "$filtered_vcf" "$snp_vcf" "$indel_vcf"

# Step 3: 统计SNP数
echo "[3/4] Counting SNPs ..."
bash /public/home/stu_tongyihan/project/7.solo_and_wild_vs_sunset/scripts/count_snp.sh \
    "$snp_vcf" "$snp_count"

# Step 4: 统计Indel指标
echo "[4/4] Calculating Indel statistics ..."
bash /public/home/stu_tongyihan/project/7.solo_and_wild_vs_sunset/scripts/indel_stats.sh \
    "$indel_vcf" "$indel_stats"

echo "=== SNP/Indel Pipeline Finished Successfully ==="
echo "Results:"
echo " - $filtered_vcf"
echo " - $snp_vcf / $indel_vcf"
echo " - $snp_count"
echo " - $indel_stats"
