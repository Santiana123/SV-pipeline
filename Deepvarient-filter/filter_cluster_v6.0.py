#!/usr/bin/env python3
import sys
import argparse
import math
import gzip
import re

def parse_args():
    parser = argparse.ArgumentParser(
        description="[v4.0] 自动统计区域密度并进行动态 Cluster 过滤 (SDR vs Autosomes)"
    )
    parser.add_argument("input_vcf", help="输入的 VCF 文件 (.vcf 或 .vcf.gz)")
    
    parser.add_argument("--sdr-chrom", type=str, default="Chr1_RagTag",
                        help="指定 SDR 所在的染色体名称 (默认: Chr1_RagTag)")
    
    parser.add_argument("-p", "--p-value", type=float, default=0.05,
                        help="统计学切断概率 (默认: 0.05)。小于此概率的间距被视为非随机簇。")
    
    parser.add_argument("--fai", type=str, default="",
                        help="基因组fai索引文件路径（优先使用此文件中的染色体长度）")
    
    # 新增参数：允许自定义阈值下限（默认10bp，也可通过命令行调整）
    parser.add_argument("--min-threshold", type=int, default=10,
                        help="过滤距离阈值的最小值 (bp)，避免阈值为0，默认10bp")
    
    return parser.parse_args()

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    else:
        return open(path, "r", encoding="utf-8", errors="replace")

def get_chrom_lengths_from_fai(fai_path):
    """从fai文件中提取染色体长度"""
    chrom_lens = {}
    if not fai_path:
        return chrom_lens
    
    try:
        with open(fai_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    chrom = parts[0]
                    try:
                        length = int(parts[1])
                        chrom_lens[chrom] = length
                    except ValueError:
                        continue
        print(f"[INFO] 从FAI文件加载了 {len(chrom_lens)} 条染色体长度信息", file=sys.stderr)
    except FileNotFoundError:
        print(f"[ERROR] FAI文件 {fai_path} 不存在！", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] 读取FAI文件失败: {str(e)}", file=sys.stderr)
        sys.exit(1)
    
    return chrom_lens

def get_chrom_lengths_from_vcf(vcf_path):
    """从 VCF Header (##contig) 中提取染色体长度"""
    chrom_lens = {}
    print(f"[Pass 1] 正在扫描 VCF Header 获取染色体长度...", file=sys.stderr)
    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("##contig"):
                # 解析 ID=xxx,length=123
                m_id = re.search(r'ID=([^,>]+)', line)
                m_len = re.search(r'length=(\d+)', line)
                if m_id and m_len:
                    chrom_lens[m_id.group(1)] = int(m_len.group(1))
            elif not line.startswith("#"):
                break # Header 结束
    return chrom_lens

def get_combined_chrom_lengths(args):
    """整合FAI和VCF Header中的染色体长度（FAI优先）"""
    # 1. 优先从FAI文件读取
    chrom_lens = get_chrom_lengths_from_fai(args.fai)
    
    # 2. 补充VCF Header中的长度（仅补充FAI中没有的）
    vcf_chrom_lens = get_chrom_lengths_from_vcf(args.input_vcf)
    for chrom, length in vcf_chrom_lens.items():
        if chrom not in chrom_lens:
            chrom_lens[chrom] = length
    
    return chrom_lens

def calculate_stats(vcf_path, chrom_lens, sdr_chrom):
    """第一遍扫描：统计各区域变异数"""
    print(f"[Pass 1] 正在统计变异数量以计算密度...", file=sys.stderr)
    
    count_sdr = 0
    count_other = 0
    
    # 如果 Header/FAI 里没找到长度，记录最大位置作为替补
    max_pos_sdr = 0
    max_pos_other = 0
    
    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"): continue
            
            fields = line.split('\t')
            if len(fields) < 2: continue
            
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                continue

            if chrom == sdr_chrom:
                count_sdr += 1
                if pos > max_pos_sdr: max_pos_sdr = pos
            else:
                count_other += 1
                if pos > max_pos_other: max_pos_other = pos
    
    # --- 计算 SDR 密度 ---
    len_sdr = chrom_lens.get(sdr_chrom, 0)
    if len_sdr == 0:
        if max_pos_sdr > 0:
            print(f"[WARN] 未找到 {sdr_chrom} 长度，使用最大位点 {max_pos_sdr} 估算。", file=sys.stderr)
            len_sdr = max_pos_sdr
        else:
            # 极端情况：VCF里SDR没有变异，避免除零
            len_sdr = 1000000 
            
    het_sdr = count_sdr / len_sdr if len_sdr > 0 else 0

    # --- 计算 Background (Other) 密度 ---
    # 计算所有非 SDR 染色体的总长度
    len_total_other = 0
    for c, l in chrom_lens.items():
        if c != sdr_chrom:
            len_total_other += l
    
    if len_total_other == 0:
        # 如果仍未找到长度，使用经验值
        print(f"[WARN] 未找到其他染色体长度，无法精确计算背景密度。", file=sys.stderr)
        print(f"[WARN] 将使用经验值 (假设其余基因组大小 300Mb) 进行估算。", file=sys.stderr)
        len_total_other = 372000000 # 番木瓜大致基因组大小
    
    het_other = count_other / len_total_other if len_total_other > 0 else 0
    
    return het_sdr, het_other

# 修改阈值计算函数：加入下限控制
def calculate_dist_threshold(het, p_value, min_threshold):
    if het <= 1e-9: 
        return min_threshold  # 原逻辑返回0，现在返回下限值
    calc_thresh = int(-math.log(1 - p_value) / het)
    # 取计算值和下限值的较大者，保证阈值≥min_threshold
    return max(calc_thresh, min_threshold)

def process_vcf(args):
    # --- 步骤 1: 获取长度与统计密度 ---
    chrom_lens = get_combined_chrom_lengths(args)
    het_sdr, het_other = calculate_stats(args.input_vcf, chrom_lens, args.sdr_chrom)
    
    # 如果背景杂合度太低（可能是空数据），给一个极小值避免报错
    if het_other == 0: het_other = 0.000001
    if het_sdr == 0: het_sdr = 0.000001

    # --- 步骤 2: 计算阈值（传入下限参数） ---
    thresh_sdr = calculate_dist_threshold(het_sdr, args.p_value, args.min_threshold)
    thresh_other = calculate_dist_threshold(het_other, args.p_value, args.min_threshold)
    
    print(f"\n[INFO] ==== 自动统计结果 ====", file=sys.stderr)
    print(f"[SDR区域] {args.sdr_chrom}:", file=sys.stderr)
    print(f"   - 观测杂合度 (Het) : {het_sdr:.6f} (约每 {int(1/het_sdr)} bp 一个变异)", file=sys.stderr)
    print(f"   - 过滤距离阈值     : {thresh_sdr} bp (P={args.p_value}, 下限={args.min_threshold}bp)", file=sys.stderr)
    
    print(f"[其他区域] Autosomes (Background):", file=sys.stderr)
    print(f"   - 观测杂合度 (Het) : {het_other:.6f} (约每 {int(1/het_other)} bp 一个变异)", file=sys.stderr)
    print(f"   - 过滤距离阈值     : {thresh_other} bp (P={args.p_value}, 下限={args.min_threshold}bp)", file=sys.stderr)
    print(f"============================\n", file=sys.stderr)

    # --- 步骤 3: 正式过滤 ---
    buffer = []
    kept_count = 0
    removed_count = 0
    
    with open_vcf(args.input_vcf) as f:
        for line in f:
            if line.startswith("#"):
                print(line, end='')
                continue
            
            fields = line.split('\t')
            try:
                chrom = fields[0]
                pos = int(fields[1])
            except ValueError:
                continue
            
            # 动态选择阈值
            current_thresh = thresh_sdr if chrom == args.sdr_chrom else thresh_other
            
            current_var = {'line': line, 'chrom': chrom, 'pos': pos, 'keep': True}
            new_buffer = []
            
            for prev_var in buffer:
                if prev_var['chrom'] != chrom:
                    if prev_var['keep']:
                        print(prev_var['line'], end='')
                        kept_count += 1
                    else:
                        removed_count += 1
                    continue # 释放旧染色体 buffer

                # 同染色体比较
                dist = pos - prev_var['pos']
                
                if dist > current_thresh:
                    # 安全
                    if prev_var['keep']:
                        print(prev_var['line'], end='')
                        kept_count += 1
                    else:
                        removed_count += 1
                else:
                    # 冲突：连坐
                    prev_var['keep'] = False
                    current_var['keep'] = False
                    new_buffer.append(prev_var)
            
            new_buffer.append(current_var)
            buffer = new_buffer
            
    # 清空 buffer
    for var in buffer:
        if var['keep']:
            print(var['line'], end='')
            kept_count += 1
        else:
            removed_count += 1
            
    print(f"[INFO] 过滤完成. 保留: {kept_count}, 移除: {removed_count}", file=sys.stderr)

if __name__ == "__main__":
    args = parse_args()
    process_vcf(args)
