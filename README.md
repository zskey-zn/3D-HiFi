## HAST : HiFi-C Accelerated Scaffolding Tool

## <span id="installation">Installation</span>

&emsp;&emsp;HAST has been tested and validated on servers running Linux.
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
