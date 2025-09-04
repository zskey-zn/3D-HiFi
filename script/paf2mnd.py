#!/usr/bin/env python3
"""
Created on Thu Aug 14 14:56:41 2025

@author: zhengshang@frasergen.com zhengshang-zn@qq.com

基本转换：
python paf2mnd.py reads.paf contacts.mnd
选取更高的identity的结果
python paf2mnd.py reads.paf contacts.mnd -m 0.85
处理超大型文件：
python paf2mnd.py huge.paf huge_contacts.mnd -m 0.75 -c 5000000 -w 16
"""

from itertools import combinations
import multiprocessing
import os
import time
from collections import defaultdict
import resource
import sys
import tempfile
import shutil
import argparse

def strand_change(record):
    """方向转换函数"""
    return 0 if record == "+" else 16

def process_read_group(group, min_identity=0.75):
    """处理单个read ID的分组数据，生成MND记录"""
    results = []
    fragments_set=set()
    dup_list=[]
    global group_dedup
    group_dedup=[]
    # 1. 检查是否为唯一比对 - 只保留只有一个比对记录的分组
    if len(group) > 1:
        # 检查是否有重复的完整read ID（包括位置信息）
        for record in group:
            if record[0] in fragments_set:
                dup_list.append(record[0])
            else:
                fragments_set.add(record[0])
        group_dedup=[i for i in group if i[0] not in dup_list]
        if len(group_dedup) == 0:
            # 去多重比对后group为空，跳过该组
            return []
    
    # 2. identity过滤
    filtered_alignments = []
    for record in group_dedup:
        try:
            matches = int(record[9])
            aln_len = int(record[10])
            if matches / aln_len > min_identity:
                filtered_alignments.append(record)
        except (ValueError, IndexError, ZeroDivisionError):
            continue
    
    # 3. 如果没有有效记录或只有一条记录，跳过该组
    if len(filtered_alignments) < 2:
        return []
    
    # 4. 生成所有有效组合
    combination_count = 0
    
    for pair in combinations(filtered_alignments, 2):
        try:
            strand1 = strand_change(pair[0][4])
            strand2 = strand_change(pair[1][4])
            mapq1 = pair[0][11]
            mapq2 = pair[1][11]
            chr1 = pair[0][5]
            chr2 = pair[1][5]
            
            # 计算位置（使用整数除法提高效率）
            start1, end1 = int(pair[0][7]), int(pair[0][8])
            start2, end2 = int(pair[1][7]), int(pair[1][8])
            pos1 = (start1 + end1) // 2
            pos2 = (start2 + end2) // 2
            
            # 提取基础read ID（去掉位置信息）
            base_read_id = pair[0][0].split(':')[0]
            
            # 构建MND记录
            results.append(
                f"{strand1}\t{chr1}\t{pos1}\t0\t{strand2}\t{chr2}\t{pos2}\t1\t{mapq1}\t-\t-\t{mapq2}\t-\t-\t{base_read_id}\t{base_read_id}\n"
            )
            
            combination_count += 1
        except (IndexError, ValueError):
            continue
    
    return results

def base_id_key(record):
    """提取基础read ID作为键"""
    full_id = record[0]
    return full_id.split(':')[0] if ':' in full_id else full_id

def process_chunk(chunk, min_identity, output_file):
    """处理一个数据块，确保同一基础read ID的记录连续"""
    # 按基础read ID分组
    groups = defaultdict(list)
    current_base_id = None
    current_group = []
    
    for record in chunk:
        base_id = base_id_key(record)
        
        # 当基础ID变化时，保存前一个组
        if current_base_id is not None and base_id != current_base_id:
            if current_group:
                groups[current_base_id] = current_group
                current_group = []
        
        current_group.append(record)
        current_base_id = base_id
    
    # 处理最后一个组
    if current_group:
        groups[current_base_id] = current_group
    
    # 处理分组并写入临时文件
    temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False)
    for group in groups.values():
        mnd_records = process_read_group(group, min_identity)
        for record in mnd_records:
            temp_file.write(record)
    temp_file.close()
    
    return temp_file.name

def memory_optimized_paf_processing(args):
    """内存优化的PAF处理流程，利用有序特性确保组完整性"""
    start_time = time.time()
    
    # 从参数对象中获取值
    paf_file = args.input
    mnd_file = args.output
    min_identity = args.min_identity
    chunk_size = args.chunk_size
    max_workers = args.max_workers
    
    # 自动确定工作进程数
    if max_workers <= 0:
        max_workers = max(1, os.cpu_count() // 2)
    
    # 使用多进程池
    pool = multiprocessing.Pool(processes=max_workers)
    
    # 分块读取和处理PAF文件
    chunk = []
    current_base_id = None
    total_lines = 0
    processed_lines = 0
    file_size = os.path.getsize(paf_file)
    processed_bytes = 0
    temp_files = []
    
    with open(paf_file, 'r') as f:
        for line in f:
            total_lines += 1
            processed_bytes += len(line.encode('utf-8'))
            
            # 解析记录
            parts = line.strip().split('\t')
            if len(parts) < 12:
                continue
                
            # 提取基础ID
            base_id = base_id_key(parts)
            
            # 当基础ID变化且块大小达到时，处理当前块
            if current_base_id is not None and base_id != current_base_id and len(chunk) >= chunk_size:
                # 提交当前块进行处理
                temp_files.append(pool.apply_async(process_chunk, (chunk.copy(), min_identity, mnd_file)))
                chunk = []
            
            chunk.append(parts)
            current_base_id = base_id
            
        # 处理最后一块
        if chunk:
            temp_files.append(pool.apply_async(process_chunk, (chunk.copy(), min_identity, mnd_file)))
    
    # 关闭进程池并等待所有任务完成
    pool.close()
    sys.stderr.write("\n等待工作进程完成...\n")
    sys.stderr.flush()
    pool.join()
    
    # 收集处理结果
    temp_file_names = [result.get() for result in temp_files]
    
    # 合并所有临时文件
    sys.stderr.write("合并临时文件...\n")
    sys.stderr.flush()
    with open(mnd_file, 'w') as out_f:
        for temp_file in temp_file_names:
            with open(temp_file, 'r') as in_f:
                shutil.copyfileobj(in_f, out_f)
            os.unlink(temp_file)
    
    # 统计结果
    total_mnd_records = 0
    with open(mnd_file, 'r') as f:
        for line in f:
            total_mnd_records += 1
    
    end_time = time.time()
    
    sys.stderr.write(f"\n处理完成: 总耗时 {end_time - start_time:.2f} 秒\n")
    sys.stderr.write(f"输入行数: {total_lines:,}\n")
    sys.stderr.write(f"生成 MND 记录数: {total_mnd_records:,}\n")
    sys.stderr.write(f"峰值内存使用: {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024:.1f} MB\n")
    sys.stderr.flush()
    
    return total_mnd_records
    
def main():
    # 创建参数解析器
    parser = argparse.ArgumentParser(
        description='将PAF比对文件转换为MND格式，用于Hi-C分析',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # 添加必需参数
    parser.add_argument('input', help='输入PAF文件路径')
    parser.add_argument('output', help='输出MND文件路径')
    
    # 添加可选参数
    parser.add_argument('-m', '--min-identity', type=float, default=0.75,
                        help='最小identity阈值 (0.0-1.0)')
    parser.add_argument('-c', '--chunk-size', type=int, default=1000000,
                        help='每个处理块的行数')
    parser.add_argument('-w', '--max-workers', type=int, default=0,
                        help='最大工作进程数 (0=自动检测)')
    
    # 解析参数
    args = parser.parse_args()
    
    # 验证参数
    if not os.path.exists(args.input):
        sys.stderr.write(f"错误: 输入文件不存在 - {args.input}\n")
        sys.exit(1)
        
    if args.min_identity < 0 or args.min_identity > 1:
        sys.stderr.write("错误: 最小比对质量必须在0.0和1.0之间\n")
        sys.exit(1)
        
    if args.chunk_size < 1000:
        sys.stderr.write("警告: 块大小过小可能导致效率低下，建议至少1000行\n")
    
    # 设置默认工作进程数
    if args.max_workers <= 0:
        args.max_workers = max(1, os.cpu_count() // 2)
    
    # 处理PAF文件
    memory_optimized_paf_processing(args)

if __name__ == "__main__":
    # 提高资源限制（针对大型文件）
    try:
        resource.setrlimit(resource.RLIMIT_NOFILE, (65536, 65536))
    except:
        pass
    
    # 设置递归深度限制（针对深层数据结构）
    sys.setrecursionlimit(10000)
    
    main()

