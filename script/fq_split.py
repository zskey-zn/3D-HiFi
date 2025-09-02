# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 14:42:30 2025
Fixed on Tue Sep  2 09:03:33 CST 2025
@author: zhengshang@frasergen.com zhengshang-zn@qq.com

# 针对大内存机器（128GB+）
python3 fq_script.py --enzyme_site GATC --fq_in huge.fq.gz --fq_out output.fq.gz \
                 --chunk_size 10000 --threads 32

# 针对普通服务器（32GB内存）
python3 fq_script.py --enzyme_site GATC --fq_in large.fq --fq_out output.fq \
                 --chunk_size 2000 --threads 16

"""
import re
import argparse
import multiprocessing
from functools import partial
import time
import os
import subprocess
import sys

def process_record(lines, compiled_pattern, min_length):
    """处理单个FASTQ记录（4行文本），切割并过滤片段"""
    if len(lines) != 4:
        return []
    
    header, seq, plus, qual = lines
    fragments = []
    header=header.split()[0]
    plus="+"
    
    # 查找酶切位点
    enzyme_site_pos = [m.start() for m in compiled_pattern.finditer(seq)]
    
    pos_start = 0
    for pos_end in enzyme_site_pos:
        frag_length = pos_end - pos_start
        if frag_length >= min_length:
            # 创建片段
            frag_seq = seq[pos_start:pos_end]
            frag_qual = qual[pos_start:pos_end]
            frag_header = f"{header}:{pos_start+1}-{pos_end}"
            
            fragments.append(f"{frag_header}\n{frag_seq}\n{plus}\n{frag_qual}\n")
        pos_start = pos_end
    
    # 处理最后一段
    frag_length = len(seq) - pos_start
    if frag_length >= min_length:
        frag_seq = seq[pos_start:]
        frag_qual = qual[pos_start:]
        frag_header = f"{header}:{pos_start+1}-{len(seq)}"
        
        fragments.append(f"{frag_header}\n{frag_seq}\n{plus}\n{frag_qual}\n")
    
    return fragments

def read_fastq_chunks(input_handle, chunk_size=1000):
    """从文件句柄中读取FASTQ记录块"""
    chunk = []
    lines = []
    
    for line in input_handle:
        lines.append(line.strip())
        if len(lines) == 4:  # 一个完整的FASTQ记录
            chunk.append(lines)
            lines = []
            
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
    
    # 处理最后可能不完整的块
    if chunk:
        yield chunk

def main():
    parser = argparse.ArgumentParser(
        description='Split FASTQ sequences at enzyme sites and filter by length',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--enzyme_site', required=True,
                        help='Enzyme recognition site (e.g., GATC)')
    parser.add_argument('--fq_in', required=True,
                        help='Input FASTQ file (gzipped if ends with .gz)')
    parser.add_argument('--fq_out', required=True,
                        help='Output FASTQ file (will gzip if ends with .gz)')
    parser.add_argument('--min_length', type=int, default=50,
                        help='Minimum fragment length to output')
    parser.add_argument('--threads', type=int, default=multiprocessing.cpu_count(),
                        help='Number of parallel threads to use')
    parser.add_argument('--chunk_size', type=int, default=1000,
                        help='Number of records per processing chunk')

    args = parser.parse_args()

    # 验证参数
    if args.min_length < 1:
        raise ValueError("Minimum length must be at least 1")
    if args.threads < 1:
        raise ValueError("Thread count must be at least 1")
    if args.chunk_size < 1:
        raise ValueError("Chunk size must be at least 1")

    print("Processing parameters:")
    print(f"  Enzyme site: {args.enzyme_site}")
    print(f"  Input file: {args.fq_in}")
    print(f"  Output file: {args.fq_out}")
    print(f"  Minimum fragment length: {args.min_length}")
    print(f"  Threads: {args.threads}")
    print(f"  Chunk size: {args.chunk_size}")
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")

    # 预编译正则表达式
    compiled_pattern = re.compile(re.escape(args.enzyme_site))
    
    # 检查是否安装了pigz/unpigz
    pigz_available = os.system("which unpigz > /dev/null 2>&1") == 0
    
    # 打开输入文件
    if args.fq_in.endswith('.gz') and pigz_available:
        print("Using unpigz for parallel decompression...")
        proc_in = subprocess.Popen(['unpigz', '-c', args.fq_in], stdout=subprocess.PIPE, text=True)
        in_handle = proc_in.stdout
    elif args.fq_in.endswith('.gz'):
        print("unpigz not found, using gzip...")
        import gzip
        in_handle = gzip.open(args.fq_in, 'rt')
    else:
        in_handle = open(args.fq_in, 'r')

    # 打开输出文件
    if args.fq_out.endswith('.gz') and pigz_available:
        print("Using pigz for parallel compression...")
        # 直接打开输出文件，稍后使用subprocess处理压缩
        out_file = open(args.fq_out, 'wb')
        proc_out = subprocess.Popen(['pigz', '-c', '-p', str(args.threads)], 
                                   stdin=subprocess.PIPE, stdout=out_file)
        out_handle = proc_out.stdin
        use_binary_output = True
    elif args.fq_out.endswith('.gz'):
        print("pigz not found, using gzip...")
        import gzip
        out_handle = gzip.open(args.fq_out, 'wt')
        use_binary_output = False
    else:
        out_handle = open(args.fq_out, 'w')
        use_binary_output = False

    # 创建进程池
    with multiprocessing.Pool(processes=args.threads) as pool:
        process_func = partial(process_record,
                              compiled_pattern=compiled_pattern,
                              min_length=args.min_length)

        # 统计变量
        total_records = 0
        total_fragments = 0
        last_report = time.time()

        # 处理记录块
        for chunk in read_fastq_chunks(in_handle, args.chunk_size):
            # 处理当前块
            results = pool.map(process_func, chunk)
            
            # 写入输出
            for fragments in results:
                for frag in fragments:
                    if use_binary_output:
                        # 对于二进制输出，需要编码为字节
                        out_handle.write(frag.encode('utf-8'))
                    else:
                        # 对于文本输出，直接写入字符串
                        out_handle.write(frag)
            
            total_records += len(chunk)
            total_fragments += sum(len(frags) for frags in results)

            # 进度报告
            current_time = time.time()
            if current_time - last_report > 60 or total_records % 10000 == 0:
                frag_per_rec = total_fragments / total_records if total_records else 0
                print(f"Processed {total_records:,} records | "
                      f"Fragments: {total_fragments:,} | "
                      f"Avg fragments/record: {frag_per_rec:.2f} | "
                      f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
                last_report = current_time

    # 清理资源
    in_handle.close()
    
    # 如果是通过subprocess打开的输入，需要等待进程结束
    if args.fq_in.endswith('.gz') and pigz_available:
        proc_in.wait()
    
    # 关闭输出句柄
    if use_binary_output:
        out_handle.close()
        proc_out.wait()
        out_file.close()
    else:
        out_handle.close()

    # 最终报告
    frag_per_rec = total_fragments / total_records if total_records else 0
    print("\nProcessing completed!")
    print(f"Total records processed: {total_records:,}")
    print(f"Total fragments generated: {total_fragments:,}")
    print(f"Average fragments per record: {frag_per_rec:.2f}")
    print(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Results saved to: {args.fq_out}")

if __name__ == "__main__":
    main()
