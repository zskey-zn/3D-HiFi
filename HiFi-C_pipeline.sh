#!/bin/bash

set -e  # 任何命令失败则退出

# 帮助信息函数
usage() {
    echo "使用方法: $0 -r PATH/contig.fa -i PATH/HiFi-C.fq.gz -p '-x map-hifi' -t 30 -c 10000 -e GATC -d PATH/3d-dna -o species"
    echo "  -r|--ref 		<contig genome> "
    echo "  -i|--fq_in 		<fastq file>  HiFi-C/Pore-C data "
    echo "  -p|--map_params 	<minimap2 align parameter > "
    echo "  -t|--threads	number of threads "
    echo "  -c|--chunk_size	Number of records per processing chunk "
    echo "  -e|--enzyme_site	Enzyme recognition site (e.g., GATC) "
    echo "  -d|--_3ddna_path	3ddna software path "
    echo "  -o|--output_prefix	output prefix "
    echo "  -h                 显示此帮助信息 "
    echo ""
    exit 0
}

# 初始化参数变量
ref=""
fq_in=""
map_params=""
threads=""
chunk_size=""
_3ddna_path=""
output_prefix=""

# 手动解析参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--ref)
            ref="$2"
            shift 2
            ;;
        -i|--fq_in)
            fq_in="$2"
            shift 2
            ;;
        -p|--map_params)
            map_params="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        -c|--chunk_size)
            chunk_size="$2"
            shift 2
            ;;
        -e|--enzyme_site)
            enzyme_site="$2"
            shift 2
            ;;     
        -d|--_3ddna_path)
            _3ddna_path="$2"
            shift 2
            ;;    
        -o|--output_prefix)
            output_prefix="$2"
            shift 2
            ;;            
        -h|--help)
            usage
            ;;
        *)
            echo "未知参数: $1"
            echo "使用 -h 查看帮助信息"
            exit 1
            ;;
    esac
done

# 检查必要参数
if [[ -z "$ref" || -z "$fq_in" || -z "$enzyme_site" || -z "$_3ddna_path" ]]; then
    echo "错误: 缺少必要参数"
    echo "使用 -h 查看帮助信息"
    exit 1
fi

# 检查是否为绝对路径
check_absolute_path() {
    local path="$1"
    local param_name="$2"
    if [[ "${path:0:1}" != "/" ]]; then
        echo "错误: 参数 $param_name 必须是绝对路径"
        exit 1
    fi
    if [[ ! -e "$path" ]]; then
        echo "错误: 路径 $path 不存在"
        exit 1
    fi
}

check_absolute_path "$ref" "--ref"
check_absolute_path "$fq_in" "--fq_in"
check_absolute_path "$_3ddna_path" "--_3ddna_path"

# 设置脚本目录和初始目录
script_dir=$(dirname "${BASH_SOURCE[0]}")
initial_dir=$(pwd)

# 创建工作目录并确保使用绝对路径
mkdir -p 00.fq_split
cd 00.fq_split
python3 ${script_dir}/script/fq_split.py --enzyme_site $enzyme_site --fq_in $fq_in --fq_out ${output_prefix}_enzyme_site_split.fq.gz --chunk_size "$chunk_size" --threads $threads

cd "$initial_dir"
mkdir -p 01.split_minimap
cd 01.split_minimap
minimap2 ${map_params} -c --secondary=no $ref ../00.fq_split/${output_prefix}_enzyme_site_split.fq.gz -t $threads > ${output_prefix}.paf

cd "$initial_dir"
mkdir -p 02.paf2mnd
cd 02.paf2mnd
python3 ${script_dir}/script/paf2mnd.py ../01.split_minimap/${output_prefix}.paf ${output_prefix}.mnd.txt -m 0.75 -c 5000000 -w $threads

# dedup
mkdir -p tmp
p=$(pwd)"/"
sort -T ./tmp -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n ${output_prefix}.mnd.txt > ${output_prefix}.mnd.sort.txt
awk -f ${script_dir}/script/dups.awk -v name=$p ${output_prefix}.mnd.sort.txt

cd "$initial_dir"
mkdir -p 03.3ddna
cd 03.3ddna
TMPDIR='./'

/usr/bin/bash ${_3ddna_path}/run-asm-pipeline.sh -r 0 --early-exit $ref ${initial_dir}/02.paf2mnd/merged_nodups.txt

#stat
cd "$initial_dir"
sed 's/:/\t/g' 01.split_minimap/${output_prefix}.paf |awk '{a[$1]=1}END{for(i in a){sum+=a[i]};print "mapping\t"sum}' > read.summary
awk '{a[$15]=1}END{for(i in a){sum+=a[i]};print "vaild map\t"sum}' 02.paf2mnd/${output_prefix}.mnd.sort.txt >> read.summary
awk '{a[$15]=1}END{for(i in a){sum+=a[i]};print "dedup\t"sum}' 02.paf2mnd/merged_nodups.txt >> read.summary
