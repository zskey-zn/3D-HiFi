## 3D-HiFi: A fast and accurate pipeline for chromosome-scale scaffolding from HiFi-C data

## <span id="Introduction">Introduction</span>
&emsp;&emsp;As an emerging high-resolution long-read chromosome conformation capture technique, HiFi-C currently lacks dedicated software tools. Researchers have to adapt pipelines designed for Pore-C, which may underestimate the true potential of HiFi-C data and negatively impact downstream analyses such as scaffolding. To address this issue, we have developed `3D-HiFi`, a tool specifically designed for efficient scaffolding using HiFi-C data. `3D-HiFi` also supports Pore-C data.

## <span id="Overview">Overview</span>
&emsp;&emsp;Given the substantial length of HiFi-C reads, we employed an analysis strategy involving in silico fragmentation at restriction enzyme cleavage sites prior to alignment.

![workflow](image/workflow.png)
## Table of contents
- [Dependencies](#Dependencies)
- [Installation](#installation)
- [Quick start](#quick_start)
- [Comparison](#comparison)
- [Output File](#output)
- [Get help](#help)
- [Citating](#Citaing)

## <span id="Dependencies">Dependencies</span>
Software:
- [python >=3.7](https://www.python.org/)
- [pigz](http://zlib.net/pigz/)
- [minimap2](https://github.com/lh3/minimap2)
- [seqkit](https://bioinf.shenwei.me/seqkit)
- [3d-dna](https://github.com/aidenlab/3d-dna)
- [parallel](https://www.gnu.org/software/parallel)
- [juicebox](https://github.com/aidenlab/Juicebox)
- [mummer](https://github.com/mummer4/mummer)

## <span id="Installation">Installation</span>

3D-HiFi has been tested and validated on servers running Linux.
```bash
# (1) Download 3D-HiFi from GitHub
$ git clone https://github.com/zskey-zn/3D-HiFi.git
# (2) Resolve dependencies
# We strongly recommend using conda to install dependencies. 
$ conda env create -f environment.yml
# Activate the HapHiC conda environment
$ conda activate 3D-HiFi # or: source /path/to/conda/bin/activate 3D-HiFi
# (3) Install 3d-dna
$ git clone https://github.com/aidenlab/3d-dna.git
#Give executable permissions
$ cd 3d-dna
$ chmod +x  *.sh */*
```

## <span id="quick_start">Quick start</span>

```bash
usage_example: 3D-HiFi -r PATH/contig.fa -i PATH/HiFi-C.fq.gz -p '-x map-hifi' -t 30 -c 10000 -e GATC -d PATH/3d-dna -o your_species

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
  -a, --polyploid      Enable polyploid mode can rescue collapsed contigs  (default: disabled)
```
## <span id="comparison">Comparison</span>
The following table summarizes a performance comparison of 3D-HiFi against other tools (wf-pore-c and Cphasing) across various biological datasets, detailing metrics such as valid reads, processing time, and memory usage.

<table style="width: 100%; table-layout: fixed; border-collapse: collapse;">
  <colgroup>
    <col style="width: 18%">
    <col style="width: 16%">
    <col style="width: 12%">
    <col style="width: 12%">
    <col style="width: 8%">
    <col style="width: 8%">
    <col style="width: 8%">
  </colgroup>
  <thead>
    <tr>
      <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Dataset</th>
      <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Software</th>
      <th style="padding: 8px; text-align: right; border-bottom: 2px solid #ddd;">Valid reads</th>
      <th style="padding: 8px; text-align: right; border-bottom: 2px solid #ddd;">Pairs number</th>
      <th style="padding: 8px; text-align: right; border-bottom: 2px solid #ddd;">Contacts/Reads</th>
      <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">Wall time</th>
      <th style="padding: 8px; text-align: left; border-bottom: 2px solid #ddd;">RAM</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td rowspan="3" style="padding: 8px; text-align: left; vertical-align: top; white-space: nowrap;"><a href="https://www.ncbi.nlm.nih.gov/sra/ERR14654081" target="_blank">Ceratitis_capitata</a> (27X) </td>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">3D-HiFi</td>
      <td style="padding: 8px; text-align: right;">1,185,394</td>
      <td style="padding: 8px; text-align: right;">11,508,130</td>
      <td style="padding: 8px; text-align: right;">8.75</td>
      <td style="padding: 8px; text-align: left;">32min</td>
      <td style="padding: 8px; text-align: left;">42G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">wf-pore-c</td>
      <td style="padding: 8px; text-align: right;">1,163,280</td>
      <td style="padding: 8px; text-align: right;">3,264,519</td>
      <td style="padding: 8px; text-align: right;">2.48</td>
      <td style="padding: 8px; text-align: left;">2.5h</td>
      <td style="padding: 8px; text-align: left;">39G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">Cphasing</td>
      <td style="padding: 8px; text-align: right;">1,015,477</td>
      <td style="padding: 8px; text-align: right;">5,888,352</td>
      <td style="padding: 8px; text-align: right;">4.48</td>
      <td style="padding: 8px; text-align: left;">48min</td>
      <td style="padding: 8px; text-align: left;">42G</td>
    </tr>
    <tr>
      <td rowspan="3" style="padding: 8px; text-align: left; vertical-align: top; white-space: nowrap;"><a href="https://www.ncbi.nlm.nih.gov/sra/ERR14654111" target="_blank">Anopheles_coluzzii</a> (56X) </td>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">3D-HiFi</td>
      <td style="padding: 8px; text-align: right;">2,282,982</td>
      <td style="padding: 8px; text-align: right;">64,712,163</td>
      <td style="padding: 8px; text-align: right;">27.35</td>
      <td style="padding: 8px; text-align: left;">1.2h</td>
      <td style="padding: 8px; text-align: left;">70G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">wf-pore-c</td>
      <td style="padding: 8px; text-align: right;">2,159,056</td>
      <td style="padding: 8px; text-align: right;">9,216,303</td>
      <td style="padding: 8px; text-align: right;">3.90</td>
      <td style="padding: 8px; text-align: left;">3.5h</td>
      <td style="padding: 8px; text-align: left;">65G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">Cphasing</td>
      <td style="padding: 8px; text-align: right;">2,258,817</td>
      <td style="padding: 8px; text-align: right;">50,310,732</td>
      <td style="padding: 8px; text-align: right;">21.27</td>
      <td style="padding: 8px; text-align: left;">1.2h</td>
      <td style="padding: 8px; text-align: left;">67G</td>
    </tr>
    <tr>
      <td rowspan="3" style="padding: 8px; text-align: left; vertical-align: top; white-space: nowrap;"><a href="https://www.ncbi.nlm.nih.gov/sra/ERR14275147" target="_blank">Homo_sapien</a> (26X) </td>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">3D-HiFi</td>
      <td style="padding: 8px; text-align: right;">9,198,589</td>
      <td style="padding: 8px; text-align: right;">720,830,448</td>
      <td style="padding: 8px; text-align: right;">78.28</td>
      <td style="padding: 8px; text-align: left;">1d5h</td>
      <td style="padding: 8px; text-align: left;">147G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">wf-pore-c</td>
      <td style="padding: 8px; text-align: right;">8,743,704</td>
      <td style="padding: 8px; text-align: right;">46,912,773</td>
      <td style="padding: 8px; text-align: right;">5.09</td>
      <td style="padding: 8px; text-align: left;">1d7h</td>
      <td style="padding: 8px; text-align: left;">61G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">Cphasing</td>
      <td style="padding: 8px; text-align: right;">9,144,776</td>
      <td style="padding: 8px; text-align: right;">364,842,740</td>
      <td style="padding: 8px; text-align: right;">39.62</td>
      <td style="padding: 8px; text-align: left;">17.9h</td>
      <td style="padding: 8px; text-align: left;">128G</td>
    </tr>
    <tr>
      <td rowspan="3" style="padding: 8px; text-align: left; vertical-align: top; white-space: nowrap;"><a href="https://www.ncbi.nlm.nih.gov/sra/SRR29580843" target="_blank">Plecia_longiforceps</a> (49X) </td>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">3D-HiFi</td>
      <td style="padding: 8px; text-align: right;">15,223,238</td>
      <td style="padding: 8px; text-align: right;">143,783,824</td>
      <td style="padding: 8px; text-align: right;">6.03</td>
      <td style="padding: 8px; text-align: left;">4.7h</td>
      <td style="padding: 8px; text-align: left;">96G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">wf-pore-c</td>
      <td style="padding: 8px; text-align: right;">13,474,187</td>
      <td style="padding: 8px; text-align: right;">35,116,293</td>
      <td style="padding: 8px; text-align: right;">1.47</td>
      <td style="padding: 8px; text-align: left;">11.2h</td>
      <td style="padding: 8px; text-align: left;">86G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">Cphasing</td>
      <td style="padding: 8px; text-align: right;">12,793,521</td>
      <td style="padding: 8px; text-align: right;">69,399,906</td>
      <td style="padding: 8px; text-align: right;">2.91</td>
      <td style="padding: 8px; text-align: left;">3.8h</td>
      <td style="padding: 8px; text-align: left;">88G</td>
    </tr>
    <tr>
      <td rowspan="3" style="padding: 8px; text-align: left; vertical-align: top; white-space: nowrap;"><a href="https://www.ncbi.nlm.nih.gov/sra/SRR28905076" target="_blank">Rosa_hybrida</a> (23X) </td>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">3D-HiFi</td>
      <td style="padding: 8px; text-align: right;">13,338,018</td>
      <td style="padding: 8px; text-align: right;">539,607,578</td>
      <td style="padding: 8px; text-align: right;">37.38</td>
      <td style="padding: 8px; text-align: left;">16.2h</td>
      <td style="padding: 8px; text-align: left;">82G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">wf-pore-c</td>
      <td style="padding: 8px; text-align: right;">10,956,145</td>
      <td style="padding: 8px; text-align: right;">37,981,725</td>
      <td style="padding: 8px; text-align: right;">2.63</td>
      <td style="padding: 8px; text-align: left;">1d2h</td>
      <td style="padding: 8px; text-align: left;">53G</td>
    </tr>
    <tr>
      <td style="padding: 8px; text-align: left; white-space: nowrap;">Cphasing</td>
      <td style="padding: 8px; text-align: right;">7,846,147</td>
      <td style="padding: 8px; text-align: right;">53,071,524</td>
      <td style="padding: 8px; text-align: right;">3.68</td>
      <td style="padding: 8px; text-align: left;">5.6h</td>
      <td style="padding: 8px; text-align: left;">71G</td>
    </tr>
  </tbody>
</table>

## <span id="output">Output Files</span>
Primary Output Files and Their Specifications

```
.
├── 01.split_minimap
│   ├── your_species.paf        # minimap2 result
│   ├── your_species.len        # contig size if you set --polyploid parameter
│   ├── contig.depth            # depth average contig if you set --polyploid parameter
│   ├── collapsed.contig.list   # collapsed contig list if you set --polyploid parameter
│   └── contig.dup.fasta        # rescued contig genome if you set --polyploid parameter
├── 02.paf2mnd
│   ├── your_species.mnd.txt
│   ├── your_species.mnd.sort.txt
│   ├── dups.txt
│   ├── merged_nodups.txt        # nodups mnd file
│   └── tmp
├── 02.02.paf2mnd_dup            # if you set --polyploid parameter
│   ├── your_species.mnd.txt
│   ├── your_species.mnd.dup.txt # new mnd file if you set --polyploid parameter
│   ├── your_species.mnd.sort.txt
│   ├── dups.txt
│   ├── merged_nodups.txt        # nodups mnd file
│   └── tmp
├── 03.3ddna
│   ├── contig.0.asm
│   ├── contig.0_asm.scaffold_track.txt
│   ├── contig.0_asm.superscaf_track.txt
│   ├── contig.0.cprops
│   ├── contig.0.assembly  # Input of Juicebox to manually correct
│   ├── contig.0.hic       # Input of Juicebox to manually correct
│   ├── contig.cprops
│   └── contig.mnd.txt
├── 03.3ddna_dup          # use this directory's result to manually correct if you set --polyploid parameter
│   ├── temp.contig.dup.0.asm_mnd.txt
│   ├── contig.dup.0_asm.scaffold_track.txt
│   ├── contig.dup.0_asm.superscaf_track.txt
│   ├── contig.dup.0.assembly  # Input of Juicebox to manually correct
│   └── contig.dup.0.hic       # Input of Juicebox to manually correct
└── read.summary          # reads mapping stat
```

## <span id="help">Get help</span>
### Help
For detailed instructions regarding chromosome ordering, orientation, and visualization, please see [automated_orient_visualization_pipeline](https://github.com/zskey-zn/3D-HiFi/tree/main/automated_orient_visualization).

Feel free to raise an issue at the [isssue page](https://github.com/zskey-zn/3D-HiFi/issues)

`Note:` Please ask questions on the issue page first. They are also helpful to other users.
### Contact
For addtional help, please send an email to zhengshang@frasergen.com

## <span id="Citaing">Citating</span>
If you use 3D-HiFi in your work,please cite:

