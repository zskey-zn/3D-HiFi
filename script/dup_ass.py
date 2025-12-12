# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 10:14:35 2025

@author: zhengshang@frasergen.com zhengshang-zn@qq.com
脚本说明：在原来聚类的结果里加入坍塌的contig信息
# 使用方法：
python script.py -a test.p_utg.0.asm -c test.p_utg.cprops -l collapsed.contig.list -o contig.dup.0.assembly
"""

import argparse

def main():
    parser = argparse.ArgumentParser(description='Process assembly files and generate duplicated contig assembly')
    parser.add_argument('-a', '--assembly', required=True, help='Input raw assembly file')
    parser.add_argument('-c', '--cprops', required=True, help='Input cprops file')
    parser.add_argument('-l', '--collapsed_list', required=True, help='Input collapsed contig list file')
    parser.add_argument('-o', '--output', required=True, help='Output assembly file')
    
    args = parser.parse_args()
    
    # Read collapsed contig list
    collapsed_ctg = []
    with open(args.collapsed_list, 'r') as f:
        for line in f:
            collapsed_ctg.append(line.strip())
    
    # Read cprops file and create dictionaries
    len_dict = {}
    raw_hic_dict = {}
    with open(args.cprops, 'r') as f:
        for line in f:
            parts = line.split()
            contig_id, hic_id, length = parts[0], parts[1], parts[2]
            len_dict[contig_id] = length
            raw_hic_dict[hic_id] = contig_id
    
    # Function to convert assembly record to (contig_id, orientation)
    def asm_change(record):
        if int(record) > 0:
            return (raw_hic_dict[record], "")
        else:
            return (raw_hic_dict[record[1:]], "-")
    
    # Create duplicated contig dictionaries
    count = 0
    dup_len_dict = {}
    dup_hic_dict = {}
    
    for hic_id, contig_id in raw_hic_dict.items():
        if contig_id in collapsed_ctg:
            count += 1
            dup_hic_dict[contig_id] = str(count)
            dup_len_dict[contig_id] = len_dict[contig_id]
            
            count += 1
            dup_hic_dict[f"{contig_id}_d2"] = str(count)
            dup_len_dict[f"{contig_id}_d2"] = len_dict[contig_id]
        else:
            count += 1
            dup_hic_dict[contig_id] = str(count)
            dup_len_dict[contig_id] = len_dict[contig_id]
    
    # Process assembly file
    raw_assembly = []
    with open(args.assembly, 'r') as f:
        for line in f:
            records = line.split()
            converted_records = [asm_change(record) for record in records]
            raw_assembly.append(converted_records)
    
    # Create duplicated assembly
    dup_assembly = []
    for line in raw_assembly:
        content = ''
        for contig_id, orientation in line:
            if contig_id in collapsed_ctg:
                content += f"{orientation}{dup_hic_dict[contig_id]} "
                dup_ctg_id = f"{contig_id}_d2"
                content += f"{orientation}{dup_hic_dict[dup_ctg_id]} "
            else:
                content += f"{orientation}{dup_hic_dict[contig_id]} "
        dup_assembly.append(content.strip())
    
    # Write output
    with open(args.output, 'w') as f:
        # Write cprops section
        for contig_id, hic_id in dup_hic_dict.items():
            f.write(f">{contig_id} {hic_id} {dup_len_dict[contig_id]}\n")
        
        # Write assembly section
        f.write('\n'.join(dup_assembly))

if __name__ == '__main__':
    main()