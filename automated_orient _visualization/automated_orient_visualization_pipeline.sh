#!/bin/bash
set -e  # Exit if any command fails

# Help information function
usage() {
    cat << EOF
Usage: $0 [options]

Required options:
  -c, --contig FILE        Input contig file path (absolute path required)
  -n, --chr_num INT        Number of chromosomes
  -d, --3ddna_path DIR     3d-dna tools path (absolute path required)

Optional options:
  -r, --related_genome FILE  Related genome file path (absolute path required)
  -e, --resolution INT       Resolution (default: 500000)
  -o, --output_prefix STR    Output prefix (default: output)
  -h, --help               Show this help message

Example:
  $0 -c /path/to/contig.fa -n 24 -d /path/to/3d-dna -r /path/to/reference.fa
EOF
    exit 0
}

# Initialize parameter variables
contig=""
chr_num=""
resolution="500000"  # Set default resolution
_3ddna_path=""
output_prefix="output"
related_genome=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--contig)
            contig="$2"
            shift 2
            ;;
        -n|--chr_num)
            chr_num="$2"
            shift 2
            ;;
        -r|--related_genome)
            related_genome="$2"
            shift 2
            ;;
        -e|--resolution)
            resolution="$2"
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
            echo "Unknown parameter: $1"
            usage
            exit 1
            ;;
    esac
done

# Check required parameters
if [[ -z "$contig" || -z "$chr_num" || -z "$_3ddna_path" ]]; then
    echo "Error: Missing required parameters"
    usage
    exit 1
fi

# Check if path is absolute
check_absolute_path() {
    local path="$1"
    local param_name="$2"
    if [[ "${path:0:1}" != "/" ]]; then
        echo "Error: Parameter $param_name must be an absolute path"
        exit 1
    fi
    if [[ ! -e "$path" ]]; then
        echo "Error: Path $path does not exist"
        exit 1
    fi
}

check_absolute_path "$contig" "--contig"
check_absolute_path "$_3ddna_path" "--_3ddna_path"

if [[ -n "$related_genome" ]]; then
    check_absolute_path "$related_genome" "--related_genome"
fi

# Set script directory and initial directory
script_dir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
initial_dir=$(pwd)
asm_prefix=$(basename "$contig" | sed 's/\.fa.*//')

# Create directory structure
mkdir -p 04.synteny/00.homo 04.synteny/01.ass2fasta 04.synteny/02.refsort 05.post_view

# Process when related genome is provided
if [[ -n "$related_genome" ]]; then
    cd 04.synteny/00.homo
    ln -sf "$related_genome" homo_chr.fa

    cd ../01.ass2fasta
    cat > step1.work.sh << EOF
#!/bin/bash
echo "Start time: \$(date)"
python3 "${script_dir}/script/ass2fasta.py" --assembly "${asm_prefix}.0.review.assembly" --ref "$contig" -n "$chr_num"
echo "End time: \$(date)"
touch ass2fasta.done
EOF
    chmod +x step1.work.sh

    cd ../02.refsort
    cat > step2.run_mummer.sh << EOF
#!/bin/bash
ref="../00.homo/homo_chr.fa"
query="../01.ass2fasta/${asm_prefix}.0.review.fasta"
prefix="$output_prefix"
echo "Start time: \$(date)"
nucmer --mum -l 100 -c 1000 -D 5 -t 16 "\$ref" "\$query" -p "\$prefix"
delta-filter "\${prefix}.delta" -1 > "\${prefix}.best.delta"
mummerplot -p "\$prefix" -f "\${prefix}.best.delta" -t postscript
show-coords -THrd "\${prefix}.best.delta" > "\${prefix}.best.delta.coord"
ps2pdf "\${prefix}.ps" "\${prefix}.pdf"
echo "End time: \$(date)"
touch mummer.done
EOF
    chmod +x step2.run_mummer.sh

    cat > step3.chrom_correspond.sh << EOF
#!/bin/bash
echo "Start time: \$(date)"
python3 "${script_dir}/script/mummer2table.py" -i "${output_prefix}.best.delta.coord"
python3 "${script_dir}/script/chrom_correspond.py" --correspond-table correspond_table --review-assembly "../01.ass2fasta/${asm_prefix}.0.review.assembly" --output-assembly "${asm_prefix}.0.correspond.assembly"
echo "End time: \$(date)"
touch chrom_correspond.done
EOF
    chmod +x step3.chrom_correspond.sh

    cd "$initial_dir/05.post_view"
    cat > step4.post_view.sh << EOF
#!/bin/bash
echo "Start time: \$(date)"
contig_fa="$contig"
correspond_assembly="../04.synteny/02.refsort/${asm_prefix}.0.correspond.assembly"
mnd_path="../02.paf2mnd/merged_nodups.txt"
resolution="$resolution"
sample_name="$output_prefix"
chr_count="$chr_num"
_3ddna_path="$_3ddna_path"
python3 "${script_dir}/script/post_review.py" -r "\$contig_fa" -a "\${correspond_assembly}" -e "\${resolution}" -m "\${mnd_path}" -c "\${chr_count}" -s "\${sample_name}" --3ddna_path "\${_3ddna_path}"
echo "End time: \$(date)"
touch post_view.done
EOF
    chmod +x step4.post_view.sh
else
    # Process when no related genome is provided
    cd 04.synteny/01.ass2fasta
    cat > step1.work.sh << EOF
#!/bin/bash
echo "Start time: \$(date)"
python3 "${script_dir}/script/ass2fasta.py" --assembly "${asm_prefix}.0.review.assembly" --ref "$contig" -n "$chr_num" --sort
echo "End time: \$(date)"
touch ass2fasta.done
EOF
    chmod +x step1.work.sh

    cd "$initial_dir/05.post_view"
    cat > step4.post_view.sh << EOF
#!/bin/bash
echo "Start time: \$(date)"
contig_fa="$contig"
correspond_assembly="../04.synteny/01.ass2fasta/${asm_prefix}.0.review_sort.assembly"
mnd_path="../02.paf2mnd/merged_nodups.txt"
resolution="$resolution"
sample_name="$output_prefix"
chr_count="$chr_num"
_3ddna_path="$_3ddna_path"
python3 "${script_dir}/script/post_review.py" -r "\$contig_fa" -a "\${correspond_assembly}" -e "\${resolution}" -m "\${mnd_path}" -c "\${chr_count}" -s "\${sample_name}" --3ddna_path "\${_3ddna_path}"
echo "End time: \$(date)"
touch post_view.done
EOF
    chmod +x step4.post_view.sh
fi

echo "Script setup completed. Please check the generated script files and run them in order."
