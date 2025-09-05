# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 10:33:50 2025

@author: zhengshang@frasergen.com zhengshang-zn@qq.com
脚本说明：处理坐标文件，提取特定染色体的信息并生成对应表
# 指定输入文件，其他使用默认
python mummer2table.py -i your.best.delta.coord

# 修改默认参数
python mummer2table.py -i your.best.delta.coord -p "Chr" -o output_correspond_table
"""

import argparse
from collections import defaultdict

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(description='处理坐标文件，提取特定染色体的信息')
    parser.add_argument('-i', '--input', default='Anopheles_coluzzii.best.delta.coord',
                        help='输入文件名（默认: Anopheles_coluzzii.best.delta.coord）')
    parser.add_argument('-p', '--prefix', default='Superscaffold',
                        help='染色体前缀（默认: Superscaffold）')
    parser.add_argument('-o', '--output', default='correspond_table',
                        help='输出文件名（默认: correspond_table）')
    
    # 解析参数
    args = parser.parse_args()
    
    # 使用参数
    coord = args.input
    chr_prefix = args.prefix
    output_file = args.output

    # 读取文件并处理数据
    data = []
    try:
        with open(coord, 'r') as f:
            for line in f:
                if chr_prefix in line:
                    parts = line.strip().split('\t')
                    if len(parts) >= 11:  # 确保有足够的列
                        # 提取第9-11列 (索引8-10)
                        data.append(parts[8:11])
    except FileNotFoundError:
        print(f"错误: 找不到输入文件 '{coord}'")
        return
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return

    if not data:
        print(f"警告: 没有找到包含前缀 '{chr_prefix}' 的行")
        return

    # 统计唯一组合的出现次数
    count_dict = defaultdict(int)
    for item in data:
        key = tuple(item)
        count_dict[key] += 1

    # 转换为列表并排序
    result_list = []
    for key, count in count_dict.items():
        # 格式: 链方向, 起始位置, 染色体, 计数
        result_list.append([key[0], key[1], key[2], count])

    # 获取所有唯一的染色体
    chromosomes = sorted(list(set(i[2] for i in result_list)))
    
    # 按染色体和计数排序
    result_list.sort(key=lambda x: (x[2], x[3]))

    # 为每个染色体选择计数最大的行
    chrom_data = {}
    for chr_id in chromosomes:
        # 获取该染色体的所有行
        chr_rows = [i for i in result_list if i[2] == chr_id]
        if chr_rows:
            # 选择计数最大的行（排序后最后一行）
            chrom_data[chr_id] = chr_rows[-1]

    # 准备最终输出
    output_data = []
    for chrom, item in chrom_data.items():
        strand = item[0]
        # 转换链的方向表示
        if strand == "1":
            strand = "+"
        elif strand == "-1":
            strand = "-"
        output_data.append([item[1], chrom, strand])

    # 按起始位置排序输出
    output_data.sort(key=lambda x: x[0])

    # 写入输出文件
    try:
        with open(output_file, 'w') as f:
            f.write("chrom\tsuperscaffold\tstrand\n")
            for item in output_data:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\n")
        print(f"处理完成，结果已保存到 '{output_file}'")
    except Exception as e:
        print(f"写入输出文件时出错: {e}")

if __name__ == "__main__":
    main()
        