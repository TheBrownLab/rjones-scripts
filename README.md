# rjones-scripts

This repository contains utility scripts written by Robert E. Jones

## Scripts

`gc_filter.py` - Plots (histogram) and filters genome contigs based in gc content. Usage: `python3 gc_filter.py [--filter FILTER] [--plot] in_file out_dir`

`rename_proteome_contigs.py` - Takes a directory of proteomes (with `.faa` extensions) as input. Outputs proteomes with renamed headers (`UniqueID_g<N>`) and a key tsv file containing the new headers and the original headers. Usage: `python3 rename_proteome_contigs.py in_dir out_dir`

`web_blast.py` - BLAST against an online database (nr, nt, swissprot, etc.) using blastn, blastp, blastx, tblastn, or tblastx. Usage: `python3 -u get_mito_genome.py -q QUERY -d DATABASE -b BLAST_TYPE`