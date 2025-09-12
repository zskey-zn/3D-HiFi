# Automated Order, orient and Visualization

This README provides details on applying automated ordering, orienting, and visualization to general genome assemblies.

## How to run
```
Usage: ./automated_orient_visualization_pipeline.sh [options]

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
  ./automated_orient__visualization -c /path/to/contig.fa -n 24 -d /path/to/3d-dna -r /path/to/reference.fa
```

## Description
### <span id="refsort">Order and orient whole chromosomes using a reference genome</span>
HAST has introduced a separate pipeline to order and orient whole chromosomes according to a reference genome.

To begin, you should get draft chromosomes genome by `.assembly` file which is juicebox manually corrected.
```bash
#set chromosome num
$ chr_num=3
#Can get contig.0.review.fasta(draft chromosomes genome) and contig.0.review.order(Genome-contig correspondence table).
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

### <span id="visualization">Visualization and Final genome</span>
<img width="998" height="850" alt="image" src="https://github.com/user-attachments/assets/b6e255aa-74ac-45ea-bca1-9f942ee159e9" />

Generate highly customizable contact maps and Final genome
```bash
$ correspond_assembly=contig.0.correspond.assembly
$ mnd_path=/path/02.paf2mnd/merged_nodups.txt
$ resolution=500000
$ sample_name=your_species
$ chr_count=3
$ _3ddna_path=/path/to/3d-dna
$ python3 /path/to/HAST/script/post_review.py -r contig.fa -a ${correspond_assembly} -e ${resolution} -m ${mnd_path} -c ${chr_count} -s${sample_name} --3ddna_path ${_3ddna_path}
```

Primary Output Files and Their Specifications
```
.
├── your_species.fasta        #Final genome
├── your_species.order        #Genome-contig correspondence table
├── your_species_chrom.fasta  #Final genome(Only Chromosome part)
├── your_species_chrom.order  #Chromosome-contig correspondence table(Only Chromosome part)
├── new_mnd_and_visualizer.sh #Script which generate new hic file 
├── your_species.hic          #new hic file
├── hic_viewer.sh             #Script which plot contact maps
└── your_species.pdf          #Contact maps
```

The blurriness of the heatmap can be mitigated by lowering the `-r` parameter in the `new_mnd_and_visualizer.sh` script to generate new `your_species.hic` file, then generate new contact maps
```bash
#Fixed (-r) parameter in new_mnd_and_visualizer.sh
$ sh new_mnd_and_visualizer.sh
#Can set new resolution(-rslu) to generate new contact maps; 
$ python3 /path/to/HAST/script/hic_viewer.py --hicfile your_species.hic --ref your_species_chrom.fasta --rslu 100000 --outpfix your_species --norm KR
$ Also can change normalization method. choices from 'KR', 'VC', 'VC_SQRT' and 'NONE'
$ python3 /path/to/HAST/script/hic_viewer.py --hicfile your_species.hic --ref your_species_chrom.fasta --rslu 100000 --outpfix your_species --norm VC
```

