#!/usr/bin/env python3
"""
脚本说明：处理染色体对应关系，重新排列组装序列
"""

import argparse
import pandas as pd
import sys

def chrom_correspond(correspond_table, review_assembly, output_assembly):
    '''处理染色体对应关系，重新排列组装序列'''
    # 读取对应表
    try:
        chrom_order = pd.read_csv(correspond_table, sep='\t')
    except Exception as e:
        print(f"错误: 无法读取对应表文件 '{correspond_table}': {e}")
        sys.exit(1)
    
    # 检查列名
    if list(chrom_order.columns) != ["chrom", "superscaffold", "strand"]:
        print("错误: 对应表的列名必须是 'chrom', 'superscaffold', 'strand'")
        sys.exit(1)

    # 处理 superscaffold 列
    chrom_order["superscaffold"] = chrom_order["superscaffold"].str.split("_").str.get(0)
    chrom_order["superscaffold"] = chrom_order["superscaffold"].str.replace("Superscaffold", '')
    chrom_order["superscaffold"] = chrom_order["superscaffold"].astype(int) - 1
    chrom_len = len(chrom_order)

    # 读取并处理组装文件
    try:
        with open(review_assembly, 'r') as infile, open(output_assembly, 'w') as outfile:
            chrom_list, unanchor_list, count = [], [], 0
            
            for line in infile:
                line = line.strip()
                if line.startswith(">"):
                    print(line, file=outfile)
                    continue
                
                if count < chrom_len:
                    chrom_list.append(line)
                    count += 1
                else:
                    unanchor_list.append(line)

            # 按照对应表重新排列染色体
            for record in chrom_order.itertuples():
                chrom_info = chrom_list[record.superscaffold]
                if record.strand == "-":
                    chrom_info = [str(-int(i)) for i in chrom_info.split()][::-1]
                    print(" ".join(chrom_info), file=outfile)
                else:
                    print(chrom_info, file=outfile)

            # 添加未锚定的序列
            for record in unanchor_list:
                print(record, file=outfile)
                
        print(f"处理完成，结果已保存到 '{output_assembly}'")
        
    except Exception as e:
        print(f"错误: 处理文件时出错: {e}")
        sys.exit(1)

def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(
        description="处理染色体对应关系，重新排列组装序列",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s --correspond-table correspond_table.txt --review-assembly sample.0.review.assembly --output-assembly sample.0.review.correspond.assembly
        """
    )
    
    # 添加选项参数
    parser.add_argument(
        "--correspond-table",
        required=True,
        help="对应表文件路径，列名必须是 'chrom', 'superscaffold', 'strand'"
    )
    
    parser.add_argument(
        "--review-assembly",
        required=True,
        help="待处理的组装文件路径"
    )
    
    parser.add_argument(
        "--output-assembly",
        required=True,
        help="输出文件路径"
    )
    
    # 解析参数
    args = parser.parse_args()
    
    # 调用主函数
    chrom_correspond(args.correspond_table, args.review_assembly, args.output_assembly)

if __name__ == '__main__':
    main()
