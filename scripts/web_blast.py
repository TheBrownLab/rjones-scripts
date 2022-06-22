#!/usr/bin/env python3

import argparse
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import multiprocessing as mp
from functools import partial


def get_args():
    parser = argparse.ArgumentParser(description='Script for performing web BLAST Search',
                                     usage="get_mito_genome.py -q QUERY -d DATABASE -b BLAST_TYPE [OPTIONS]",
                                     formatter_class=argparse.RawTextHelpFormatter)
    optional = parser._action_groups.pop()  # Edited this line
    required = parser.add_argument_group('required arguments')

    # Required Arguments
    required.add_argument('-q', '--query', type=str, required=True,
                          help='path to file to be used as query for BLAST search')
    required.add_argument('-d', '--database', type=str, required=True,
                          help='online database to be used. \n'
                               '(e.g. nr, nt, swissprot, etc.)')
    required.add_argument('-b', '--blast_type', type=str, required=True,
                          help='blastn, blastp, blastx, tblastn, or tblastx')

    # Optional Arguments
    optional.add_argument('-t', '--threads', type=str, default=1,
                          help='Number of threads. \n'
                               'Default=1')
    optional.add_argument('-o', '--output', type=str, default='out.fas',
                          help='Name of output file in FASTA format. \n'
                               'Default=out.fas')
    parser._action_groups.append(optional)

    return parser.parse_args()


def do_blast(seq, blast_type, db):
    blast_results = NCBIWWW.qblast(program=blast_type, database=db, sequence=seq)
    return blast_results
    # with open("my_blast.xml", "a") as save_to:
    #     save_to.write(blast_results.read())
    #     blast_results.close()
    #
    # with open('my_blast.xml', 'r') as blast_results:
    #     blast_records = NCBIXML.parse(blast_results)
    #
    #     for blast_record in blast_records:
    #         for i, description in enumerate(blast_record.alignments):
    #             print(description)
    #             # for hsp in alignment.hsps:
    #             #     print(hsp)


if __name__ == '__main__':
    open('my_blast.xml', 'w').close()
    args = get_args()
    infile = args.query
    records = SeqIO.parse(infile, format='fasta')

    with mp.Pool(processes=int(args.threads)) as p:
        # Iterable to be used in p.map. This iterates through each sequence in the input fasta file.
        _iter = [record.format(format='fasta') for record in records]

        # first argument is function to be run in parallel. 2nd and 3rd arguments are function arguments that are the
        # same for all iterations of the function
        func = partial(do_blast, blast_type=args.blast_type, db=args.database)
        results = p.map(func, iterable=_iter)

    for result in results:
        records = NCBIXML.parse(result)
        for record in records:
            for i, description in enumerate(record.alignments):
                print(description)
