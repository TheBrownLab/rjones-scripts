# rjones-scripts

This repository contains utility scripts written by Robert E. Jones

## Scripts

`rename_proteome_contigs.py` - Takes a directory of proteomes (with `.faa` extensions) as input. Outputs proteomes with renamed headers (`UniqueID_g<N>`) and a key tsv file containing the new headers and the original headers.
- Usage: `python3 rename_proteome_contigs.py in_dir out_dir`