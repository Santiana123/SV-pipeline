#!/usr/bin/env python3
import sys
import gzip
import numpy as np

if len(sys.argv) != 2:
    sys.stderr.write(f"Usage: {sys.argv[0]} <vcf_file>\n")
    sys.exit(1)

vcf_file = sys.argv[1]
dp_values = []

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    else:
        return open(path, "r")

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
        dp_found = False

        # --- 策略 1: 优先读取标准 DP (适用于 PBSV, SVIM) ---
        if "DP" in fmt:
            try:
                idx = fmt.index("DP")
                val = sample[idx]
                if val != "." and int(val) > 0:
                    dp_values.append(int(val))
                    dp_found = True
            except: pass

        # --- 策略 2: 针对 cuteSV/Sniffles2 采用带“比例平衡”的 DR + DV 统计 ---
        if not dp_found and "DR" in fmt and "DV" in fmt:
            try:
                dr_idx = fmt.index("DR")
                dv_idx = fmt.index("DV")
                dr = int(sample[dr_idx]) if sample[dr_idx] != "." else 0
                dv = int(sample[dv_idx]) if sample[dv_idx] != "." else 0
                
                # 排除纯合缺失或零深度行
                if dr + dv > 0:
                    # 引入 3 倍差值过滤逻辑 (针对 0/1 类型变异设计)
                    # 只有当两者都不为 0 且比例在 1/3 到 3 之间时才纳入统计
                    if dr > 0 and dv > 0:
                        ratio = max(dr, dv) / min(dr, dv)
                        if ratio <= 3.0:
                            dp_values.append(dr + dv)
                            dp_found = True
                    # 如果你以后想统计 1/1 类型，可以放开对 dr=0 的限制
            except: pass

if len(dp_values) == 0:
    sys.stderr.write("Error: No valid balanced depth/DP fields found.\n")
    sys.exit(1)

dp = np.array(dp_values)
median_dp = np.median(dp)
mean_dp = np.mean(dp)

print("\n==== SV Physical Depth Statistics v6.0 (Balanced-Mode) ====")
print(f"Total variants used for stats: {dp.size}")
print(f"Median Depth (P50)           : {median_dp:.2f}")
print(f"Mean Depth                   : {mean_dp:.2f}")
print(f"Min Value                    : {dp.min()}")
print(f"Max Value                    : {dp.max()}")
print("-" * 55)
print(f"Suggested minDP (0.5x Median): {max(5, int(median_dp * 0.5))}")
print(f"Suggested maxDP (2.0x Median): {int(median_dp * 2.0)}")
print("-" * 55)
print("Note: Statistics excluded variants with >3x difference between DR and DV.")
