# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:09:08 2026

@author: zhengshang@frasergen.com zhengshang-zn@qq.com

depth plot and collapsed contig report

"""
import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks

# ================= 命令行参数解析 =================
parser = argparse.ArgumentParser(description='基于深度分布鉴定塌陷contig')
parser.add_argument('depth_file', help='深度文件 (格式: contig start end depth)')
parser.add_argument('ctg_sizes', help='contig长度文件 (格式: contig length)')
parser.add_argument('-o', '--output', default='collapsed.contig.list',
                    help='输出文件名 (默认: collapsed.contig.list)')
parser.add_argument('-p', '--polyploid', type=int, default=2,
                    help='倍性，用于生成阈值等级 (默认: 2)')
parser.add_argument('-r', '--min_high_depth_ratio', type=float, default=0.5,
                    help='高深度窗口占比阈值 (默认: 0.5)')
args = parser.parse_args()

# 将命令行参数赋给变量
depth_file = args.depth_file
ctg_sizes = args.ctg_sizes
output = args.output
polyploid = args.polyploid
min_high_depth_ratio = args.min_high_depth_ratio

# ================= 其他可调参数（仍保留为脚本内常量） =================
bw_method = 0.18
peak_distance = 30
prominence_factor = 0.02
bins = 60
tail_percentile = 99.5
# =============================================

def load_contig_lengths(filepath):
    """读取 contig 长度文件，返回名称 -> 长度的字典（int）"""
    ctg_len = {}
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    name = parts[0]
                    try:
                        length = int(parts[1])
                        ctg_len[name] = length
                    except ValueError:
                        continue
    except FileNotFoundError:
        print(f"错误：找不到 contig 长度文件 {filepath}", file=sys.stderr)
        sys.exit(1)
    return ctg_len

# 1. 读取深度数据
df = pd.read_csv(depth_file, sep='\t', header=None,
                 names=['contig', 'start', 'end', 'depth'])
depth_data = df['depth'].values

# 2. 动态 x 轴上限
x_upper_limit = np.percentile(depth_data, tail_percentile)
x_range = (0, x_upper_limit + 15)

# 3. KDE 估计与寻峰
x_fit = np.linspace(x_range[0], x_range[1], 1000)
kde = gaussian_kde(depth_data, bw_method=bw_method)
pdf_fit = kde(x_fit)

peaks, _ = find_peaks(pdf_fit, distance=peak_distance,
                      prominence=pdf_fit.max() * prominence_factor)

valid_peaks = [(x_fit[i], pdf_fit[i]) for i in peaks if x_fit[i] > 2]
valid_peaks.sort(key=lambda t: t[1], reverse=True)

main_peak_depth = round(valid_peaks[0][0], 2) if valid_peaks else None
secondary_peak_depth = round(valid_peaks[1][0], 2) if len(valid_peaks) > 1 else None

print("\n" + "=" * 40)
print(f"★ 主峰深度: {main_peak_depth}")
print(f"★ 次峰深度: {secondary_peak_depth}")
print("=" * 40 + "\n")

# 4. 绘图
fig, ax = plt.subplots(figsize=(7, 5.5))
ax.hist(depth_data, bins=bins, range=x_range,
        density=True, alpha=0.4, color='skyblue', edgecolor='black')
ax.plot(x_fit, pdf_fit, color='#178a8f', linewidth=2.5)

y_max = pdf_fit.max() * 1.15
peak_colors = ['#cc1111', '#333a42']
for idx, (px, py) in enumerate(valid_peaks[:2]):
    if px > x_upper_limit:
        continue
    color = peak_colors[idx]
    v_max = min(py * (1.05 if idx == 0 else 1.25), y_max * 0.95)
    ax.vlines(x=px, ymin=0, ymax=v_max, colors=color,
              linestyles='--', linewidth=1.2)
    ax.text(px + 0.8, v_max, f'{int(round(px))}',
            color=color, fontsize=12, fontweight='bold', va='bottom')

ax.set_xlim(-2, x_upper_limit + 10)
ax.set_ylim(0, y_max)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)
ax.set_xlabel('Coverage', fontsize=14, fontweight='bold', labelpad=10)
ax.set_ylabel('Density', fontsize=14, fontweight='bold', labelpad=10)
ax.tick_params(axis='both', labelsize=12, width=1.5)
ax.grid(alpha=0.15, linestyle='--')
plt.tight_layout()
plt.savefig('depth_fit_curve_fixed.png', dpi=300)
plt.show()

# 5. 塌陷 contig 筛选
if main_peak_depth is None:
    print("未检测到主峰，无法进行塌陷筛选。")
    sys.exit(0)

print(f"正在使用主峰深度 {main_peak_depth} 进行后续过滤...")
ctg_len_dict = load_contig_lengths(ctg_sizes)

df_filter = df[df['depth'] >= 1.5 * main_peak_depth].copy()
df_filter['window_len'] = df_filter['end'] - df_filter['start']

high_depth_len = df_filter.groupby('contig').agg(
    high_depth_length=('window_len', 'sum'),
    avg_depth=('depth', 'mean')
).reset_index()

high_depth_len['full_length'] = high_depth_len['contig'].map(ctg_len_dict).astype(float)
high_depth_len.dropna(subset=['full_length'], inplace=True)
high_depth_len['full_length'] = high_depth_len['full_length'].astype(int)
high_depth_len['high_depth_ratio'] = high_depth_len['high_depth_length'] / high_depth_len['full_length']

high_depth_len_ge_50 = high_depth_len[high_depth_len['high_depth_ratio'] >= min_high_depth_ratio].copy()
print(f"高深度窗口占比 >= {min_high_depth_ratio} 的 contig 数量: {len(high_depth_len_ge_50)}")

thresholds = [(i + 0.5) * main_peak_depth for i in range(1, polyploid)]

def assign_level(avg_depth):
    for i in range(len(thresholds)-1, -1, -1):
        if avg_depth >= thresholds[i]:
            return i + 2
    return None

high_depth_len_ge_50['level'] = high_depth_len_ge_50['avg_depth'].apply(assign_level)
collapsed_contigs = high_depth_len_ge_50.dropna(subset=['level'])
collapsed_contigs['level'] = collapsed_contigs['level'].astype(int)

with open(output, 'w') as f:
    for _, row in collapsed_contigs.iterrows():
        levels_labels = [f'd{i}' for i in range(2, row['level']+1)]
        levels_str = '\t'.join(levels_labels)
        f.write(f"{row['contig']}\t{levels_str}\n")

print(f"塌陷 contig 结果已保存至 {output}，共 {len(collapsed_contigs)} 条。")
