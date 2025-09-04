## HAST : HiFi-C Accelerated Scaffolding Tool

## <span id="Introduction">Introduction</span>
&emsp;&emsp;As an emerging high-resolution long-read chromosome conformation capture technique, HiFi-C currently lacks dedicated software tools. Researchers have to adapt pipelines designed for Pore-C, which may underestimate the true potential of HiFi-C data and negatively impact downstream analyses such as scaffolding. To address this issue, we have developed `HAST`, a tool specifically designed for efficient scaffolding using HiFi-C data. `HAST` also supports Pore-C data.

## <span id="Dependencies">Dependencies</span>
Software:
- [python >=3.7](https://www.python.org/)
- [pigz](http://zlib.net/pigz/)
- [minimap2](https://github.com/lh3/minimap2)
- [3d-dna](https://github.com/aidenlab/3d-dna)
- [parallel](https://www.gnu.org/software/parallel)

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

## <span id="help">Get help</span>
### Help
Feel free to raise an issue at the [isssue page](https://github.com/zskey-zn/HAST/issues)
`Note:` Please ask questions on the issue page first. They are also helpful to other users.
### Contact
For addtional help, please send an email to zhengshang@frasergen.com

## <span id="Citaing">Citating</span>
If you use HAST in your work,please cite:

