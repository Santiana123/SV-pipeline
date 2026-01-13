#!/usr/bin/env python3
import sys
import gzip
import numpy as np

if len(sys.argv) != 2:
    sys.stderr.write(f"Usage: {sys.argv[0]} <input.vcf.gz>\n")
    sys.exit(1)

vcf_file = sys.argv[1]
dp_values = []

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r")

# --- 读取 VCF ---
print(f"Reading VCF file: {vcf_file} ...", file=sys.stderr)
with open_vcf(vcf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")
        if len(fields) < 10:
            continue

        fmt = fields[8].split(":")
        sample = fields[9].split(":")

        if "DP" not in fmt:
            continue

        try:
            dp_idx = fmt.index("DP")
            dp_str = sample[dp_idx]
            if dp_str == ".":
                continue
            dp = int(dp_str)
            dp_values.append(dp)
        except (ValueError, IndexError):
            continue

if len(dp_values) == 0:
    sys.stderr.write("No DP values found in VCF.\n")
    sys.exit(1)

dp = np.array(dp_values)

# --- 统计基础信息 ---
median_dp = np.median(dp)
mean_dp = np.mean(dp)
max_dp_actual = np.max(dp)

print("\n==== DP Statistics ====")
print(f"Total variants : {dp.size}")
print(f"Min DP         : {dp.min()}")
print(f"Mean DP        : {mean_dp:.2f}")
print(f"Max DP         : {max_dp_actual}")

# --- 详细分位数 (Quantiles) ---
quantiles_dict = {
    "P1": 0.01,
    "P5": 0.05,
    "P10": 0.10,
    "P25": 0.25,
    "Median(P50)": 0.50,
    "P75": 0.75,
    "P90": 0.90,
    "P95": 0.95,
    "P99": 0.99,
    "P99.5": 0.995
}

print("\n---- Quantiles Distribution ----")
for name, q in quantiles_dict.items():
    print(f"{name:12s}: {np.quantile(dp, q):.2f}")

# --- 阈值计算 (严格模式: 0.5x - 2.0x) ---
min_dp_strict = int(median_dp * 0.5)
max_dp_strict = int(median_dp * 2.0)

# 设定绝对下限 (防止测序深度极低时阈值变成0或1)
if min_dp_strict < 5:
    min_dp_strict = 5

print("\n==== Recommended Filters (Strict / High Confidence) ====")
print(f"Anchor (Median): {median_dp:.2f}")
print(f"Strategy       : [ 0.5 * Median, 2.0 * Median ]")
print("-" * 40)
print(f"Suggested minDP : {min_dp_strict}")
print(f"Suggested maxDP : {max_dp_strict}")
print("-" * 40)

# --- ASCII 直方图 ---
print(f"\n==== DP Histogram (Cutoff Preview: {min_dp_strict}-{max_dp_strict}) ====")
# 智能调整直方图范围：如果有极长尾，只显示到 P99 或 max_dp_strict 的 1.5 倍，避免图被拉得太扁
p99 = np.quantile(dp, 0.99)
display_max = max(max_dp_strict * 1.5, p99) 
# 如果最大值实在太大，强制截断以便观察主峰
if display_max > median_dp * 5:
    display_max = median_dp * 5

hist, bin_edges = np.histogram(dp, bins=60, range=(0, display_max))
max_hist = max(hist)
scale = 50 / max_hist if max_hist > 0 else 1

for i in range(len(hist)):
    bar_len = int(hist[i] * scale)
    bar = "*" * bar_len
    low_edge = int(bin_edges[i])
    high_edge = int(bin_edges[i+1])
    
    # 标记位置
    mark = ""
    # 标记 Median
    if low_edge <= median_dp < high_edge:
        mark += " <--- Median"
    
    # 标记 阈值边界
    if low_edge <= min_dp_strict < high_edge:
        mark += " [Min Cutoff]"
    if low_edge <= max_dp_strict < high_edge:
        mark += " [Max Cutoff]"
        
    print(f"{low_edge:3d}-{high_edge:3d} | {bar}{mark}")

print("\n==== Example bcftools command ====")
print(f"bcftools filter -i 'FORMAT/DP>={min_dp_strict} && FORMAT/DP<={max_dp_strict}' \\")
print(f"    {vcf_file} -Oz -o {vcf_file.replace('.vcf.gz', '')}.filtered.vcf.gz")
