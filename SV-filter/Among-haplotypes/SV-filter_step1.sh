#!/bin/bash

# --- 脚本路径配置 ---
# 请确保这两个脚本在对应的路径下
STATS_PY="/public/home/yuejingjing/tyh/protocol/SV-filter/SV_depth_stats_v6.0.py"
FILTER_PY="/public/home/yuejingjing/tyh/protocol/SV-filter/filter_sv_vAF.py"

# --- 阈值提取函数 ---
# 自动运行 v6.0 统计脚本并捕获建议的 min/max
get_thresholds() {
    local vcf=$1
    local stats=$(python $STATS_PY "$vcf" 2>/dev/null | grep "Suggested")
    min=$(echo "$stats" | grep "minDP" | awk -F': ' '{print $2}')
    max=$(echo "$stats" | grep "maxDP" | awk -F': ' '{print $2}')
    echo "$min $max"
}

# --- 批量处理流程 ---
tools=("cutesv" "sniffles2" "svim" "pbsv")

echo "===================================================="
echo "开始 SV 自动化过滤流程"
echo "标准: GT=0/1, VAF=0.35-0.65, 动态深度"
echo "===================================================="

for tool in "${tools[@]}"; do
    vcf="${tool}.vcf"
    if [ -f "$vcf" ]; then
        echo "[处理中] $tool ..."
        
        # 1. 动态获取上一步脚本提供的阈值
        read t_min t_max <<< $(get_thresholds "$vcf")
        
        if [ -z "$t_min" ] || [ -z "$t_max" ]; then
            echo "  - 错误: 无法从统计脚本获取 $tool 的阈值。"
            continue
        fi
        
        echo "  - 动态阈值: minDP=$t_min, maxDP=$t_max"
        
        # 2. 执行精准过滤
        python $FILTER_PY "$vcf" "$tool" "$t_min" "$t_max" > "${tool}.f.vcf"
        
        # 统计过滤前后的变异数量
        before=$(grep -v "^#" "$vcf" | wc -l)
        after=$(grep -v "^#" "${tool}.f.vcf" | wc -l)
        echo "  - 完成! 数量变化: $before -> $after"
    else
        echo "[跳过] 找不到文件: $vcf"
    fi
done

echo "===================================================="
echo "所有过滤任务已完成。输出文件为 *.f.vcf"
