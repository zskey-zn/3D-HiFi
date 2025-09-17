#!/bin/bash

set -e  # 任何命令失败则退出

# 帮助信息函数
usage() {
    echo "Usage: $0 -r PATH/contig.fa -i PATH/HiFi-C.fq.gz -p '-x map-hifi' -t 30 -c 10000 -e GATC -d PATH/3d-dna -o species"
    echo "  -h|--help           show this help message and exit"
    echo "  -r|--ref            <contig genome> "
    echo "  -i|--fq_in          <fastq file>  HiFi-C/Pore-C data "
    echo "  -p|--map_params     <minimap2 align parameter > "
    echo "  -t|--threads        Number of threads "
    echo "  -c|--chunk_size     Number of records per processing chunk. If the dataset is large, you can increase the chunk_size parameter"
    echo "  -e|--enzyme_site    Enzyme recognition site : GATC (MboI/DpnII), AAGCTT (HindIII), CATG(NlaIII)"
    echo "  -d|--_3ddna_path    3ddna software path "
    echo "  -o|--output_prefix	output prefix "
    echo "  -a|--polyploid      Enable polyploid mode can rescue collapsed contigs (default: disabled)"
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
polyploid=false  # 默认不启用多倍体模式

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
        -a|--polyploid)
            polyploid=true
            shift
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
asm_prefix=$(basename "$ref" | sed 's/\.fa.*//')
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
python3 ${script_dir}/script/paf2mnd.py ../01.split_minimap/${output_prefix}.paf ${output_prefix}.mnd.txt -m 0.75 -c 500000 -w $threads

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

cd "$initial_dir"

if [ "$polyploid" = true ]; then
    cd 01.split_minimap
    seqkit fx2tab -l -n -i $ref  > ${output_prefix}.len
    awk 'NR==FNR{contig_l[$1]=$2;genome_size+=$2}NR>FNR{contig_map[$6]+=($9-$8+1);genome_map+=($9-$8+1)}END{for(i in contig_map){print i,contig_map[i]/contig_l[i],genome_map/genome_size}}' ${output_prefix}.len ${output_prefix}.paf > contig.depth
    awk '$2>=$3*1.5{print $1}' contig.depth > collapsed.contig.list
    seqkit grep -f collapsed.contig.list $ref | seqkit seq -w0 | awk 'NR%2==1{print $1"_d2"}NR%2==0{print}' | seqkit seq $ref - | seqkit sort > contig.dup.fasta
    
    mkdir -p ../02.paf2mnd_dup
    cd ../02.paf2mnd_dup
    ln -s ../02.paf2mnd/${output_prefix}.mnd.txt .
    awk -v OFS="\t" 'NR==FNR{a[$1]=1}NR>FNR&&a[$2]&&a[$6]{$2=$2"_d2";$6=$6"_d2";print $0}' ../01.split_minimap/collapsed.contig.list ${output_prefix}.mnd.txt > ${output_prefix}.mnd.dup.txt
    awk -v OFS="\t" 'NR==FNR{a[$1]=1}NR>FNR&&a[$2]{$2=$2"_d2";print $0}' ../01.split_minimap/collapsed.contig.list ${output_prefix}.mnd.txt >> ${output_prefix}.mnd.dup.txt
    awk -v OFS="\t" 'NR==FNR{a[$1]=1}NR>FNR&&a[$6]{$6=$6"_d2";print $0}' ../01.split_minimap/collapsed.contig.list ${output_prefix}.mnd.txt >> ${output_prefix}.mnd.dup.txt
    cat ${output_prefix}.mnd.txt >> ${output_prefix}.mnd.dup.txt
    # dedup
    mkdir -p tmp
    p=$(pwd)"/"
    sort -T ./tmp -m -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n ${output_prefix}.mnd.dup.txt > ${output_prefix}.mnd.sort.txt
    awk -f ${script_dir}/script/dups.awk -v name=$p ${output_prefix}.mnd.sort.txt
    
    mkdir -p ../03.3ddna_dup
    cd ../03.3ddna_dup
    python3 ${script_dir}/script/dup_ass.py --assembly ../03.3ddna/${asm_prefix}.0.asm --cprops ../03.3ddna/${asm_prefix}.0.cprops --collapsed_list ../01.split_minimap/collapsed.contig.list --output contig.dup.0.assembly
    /usr/bin/bash ${_3ddna_path}/visualize/run-assembly-visualizer.sh -p true contig.dup.0.assembly ../02.paf2mnd_dup/merged_nodups.txt

fi

#stat
cd "$initial_dir"
sed 's/:/\t/g' 01.split_minimap/${output_prefix}.paf |awk '{a[$1]=1}END{for(i in a){sum+=a[i]};print "mapping\t"sum}' > read.summary
awk '{a[$15]=1}END{for(i in a){sum+=a[i]};print "vaild map\t"sum}' 02.paf2mnd/${output_prefix}.mnd.sort.txt >> read.summary
awk '{a[$15]=1}END{for(i in a){sum+=a[i]};print "dedup\t"sum}' 02.paf2mnd/merged_nodups.txt >> read.summary
