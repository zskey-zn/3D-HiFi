## HAST : HiFi-C Accelerated Scaffolding Tool

## <span id="Introduction">Introduction</span>
&emsp;&emsp;As an emerging high-resolution long-read chromosome conformation capture technique, HiFi-C currently lacks dedicated software tools. Researchers have to adapt pipelines designed for Pore-C, which may underestimate the true potential of HiFi-C data and negatively impact downstream analyses such as scaffolding. To address this issue, we have developed `HAST`, a tool specifically designed for efficient scaffolding using HiFi-C data. `HAST` also supports Pore-C data.

## <span id="Overview">Overview</span>
&emsp;&emsp;Given the substantial length of HiFi-C reads, we employed an analysis strategy involving in silico fragmentation at restriction enzyme cleavage sites prior to alignment.

![workflow](image/workflow.png)
## Table of contents
- [Dependencies](#Dependencies)
- [Installation](#installation)
- [Quick start](#quick_start)
- [Output File](#output)
- [Order and orient whole chromosomes using a reference genome](#refsort)
- [Get help](#help)
- [Citating](#Citaing)
## <span id="Dependencies">Dependencies</span>
Software:
- [python >=3.7](https://www.python.org/)
- [pigz](http://zlib.net/pigz/)
- [minimap2](https://github.com/lh3/minimap2)
- [3d-dna](https://github.com/aidenlab/3d-dna)
- [parallel](https://www.gnu.org/software/parallel)
- [juicebox](https://github.com/aidenlab/Juicebox)
- [mummer](https://github.com/mummer4/mummer)

## <span id="Installation">Installation</span>

HAST has been tested and validated on servers running Linux.
```bash
# (1) Download HAST from GitHub
$ git clone https://github.com/zskey-zn/HAST.git
# (2) Resolve dependencies
# We strongly recommend using conda to install dependencies. 
$ conda env create -f environment.yml
# Activate the HapHiC conda environment
$ conda activate hast # or: source /path/to/conda/bin/activate hast
# (3) Install 3d-dna
$ git clone https://github.com/aidenlab/3d-dna.git
#Give executable permissions
$ cd 3d-dna
$ chmod +x  *.sh */*
```

## <span id="quick_start">Quick start</span>

```bash
usage_example: HiFi-C_pipeline.sh -r PATH/contig.fa -i PATH/HiFi-C.fq.gz -p '-x map-hifi' -t 30 -c 10000 -e GATC -d PATH/3d-dna -o your_species

options:
  -h, --help           show this help message and exit
  -r, --ref            contig genome
  -i, --fq_in          <fastq file>  HiFi-C/Pore-C data
  -p, --map_params     <minimap2 align parameter> (if your data is Pore-C,set "-x map-ont" )
  -t, --threads        number of threads
  -c, --chunk_size     Number of records per processing chunk, If the dataset is large, you can increase the `chunk_size` parameter.
  -e, --enzyme_site    Enzyme recognition site: GATC (MboI/DpnII), AAGCTT (HindIII), CATG(NlaIII)
  -d, --_3ddna_path    3ddna software path
  -o, --output_prefix  output prefix
```

## <span id="output">Output Files</span>
Primary Output Files and Their Specifications

```
.
├── 02.paf2mnd
│   ├── your_species.mnd.sort.txt
│   ├── your_species.mnd.txt
│   ├── dups.txt
│   ├── merged_nodups.txt  #nodups mnd file
│   └── tmp
├── 03.3ddna
│   ├── contig.0.asm
│   ├── contig.0_asm.scaffold_track.txt
│   ├── contig.0_asm.superscaf_track.txt
│   ├── contig.0.cprops
│   ├── contig.0.assembly  #Input of Juicebox to manually correct
│   ├── contig.0.hic       #Input of Juicebox to manually correct
│   ├── contig.cprops
│   └── contig.mnd.txt  
└── read.summary          #reads mapping stat
```
## <span id="refsort">Order and orient whole chromosomes using a reference genome</span>
HAST has introduced a separate pipeline to order and orient whole chromosomes according to a reference genome.

To begin, you should get draft chromosomes genome by `.assembly` file which is juicebox manually corrected.
```bash
#set chromosome num
$ chr_num=3
#Can get contig.0.review.fasta(draft chromosomes genome) and contig.0.review.order(Chromosome-contig correspondence table).
$ python3 /path/to/HAST/script/ass2fasta.py --assembly contig.0.review.assembly --ref contig.fa -n $chr_num
#If don't have closely related chromosomes genome, can sort in descending order by Chromosome length by adding "--sort" parameter，will obtain the contig.0.review_sort.assembly file(sorted assembly file) as an additional output.
$ python3 /path/to/HAST/script/ass2fasta.py --assembly contig.0.review.assembly --ref contig.fa -n $chr_num --sort
```

Then use [mummer](https://github.com/mummer4/mummer) align draft chromosomes genome to  a chromosome-level reference genome. The reference genome can be from the same species or a closely related one:
```bash
$ ref=/path/to/closely_related.genome.fa
$ query=contig.0.review.fasta
$ prefix=your_species
$ nucmer --mum -l 100 -c 1000 -D 5 -t 16 $ref $query -p $prefix  # *** you should change the paramenter when your genome is too large ,like  -l 500
$ delta-filter  ${prefix}.delta -1 > ${prefix}.best.delta  # if has noisy ,you can add -i -o -l paramenter ,like  -i 95 -o 95 -l 150
$ mummerplot -p $prefix -f ${prefix}.best.delta -t postscript
$ show-coords -THrd ${prefix}.best.delta > ${prefix}.best.delta.coord
$ /usr/bin/ps2pdf ${prefix}.ps ${prefix}.pdf
#Generate the correspondence table(correspond_table), can fixed by ${prefix}.pdf
$ python3 /path/to/HAST/script/mummer2table.py -i ${prefix}.best.delta.coord
```

Finally, generate the final file based on the correspondence table.
```bash
$ python3 /path/to/HAST/script/chrom_correspond.py --correspond-table correspond_table --review-assembly contig.0.review.assembly --output-assembly contig.0.correspond.assembly
```



## <span id="help">Get help</span>
### Help
Feel free to raise an issue at the [isssue page](https://github.com/zskey-zn/HAST/issues)

`Note:` Please ask questions on the issue page first. They are also helpful to other users.
### Contact
For addtional help, please send an email to zhengshang@frasergen.com

## <span id="Citaing">Citating</span>
If you use HAST in your work,please cite:

