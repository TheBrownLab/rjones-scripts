#! /usr/bin/env python3
'''
Replicates missing data in simulated matrix.
'''

import argparse as ap
from Bio import SeqIO

def insert_missing_data(real, sim, outfile):
    '''
    Replicates missing data in simulated matrix.

    :param real: path to input real matrix
    :type real: str
    :param sim: path to input simulated matrix
    :type sim: str
    :param outfile: path to output file
    :type outfile: str
    '''
    # Parse input files
    real_records = [record for record in SeqIO.parse(real, 'fasta')]
    sim_records = [record for record in SeqIO.parse(sim, 'phylip-relaxed')]
    # Sort based on taxa name
    real_records = sorted(real_records, key=lambda x: x.id)
    sim_records = sorted(sim_records, key=lambda x: x.id)
    # Make sure same number of taxa in each
    assert len(real_records) == len(sim_records)

    # Initialize list of records
    out_records = []

    # Iterate through SeqRecords
    for real_record, sim_record in zip(real_records, sim_records):
        assert real_record.id == sim_record.id
        assert len(real_record.seq) == len(sim_record.seq)
        for i, base in enumerate(real_record.seq):
            if base == '-':
                sim_record.seq = sim_record.seq[:i] + '-' + sim_record.seq[i+1:]
        out_records.append(sim_record)

    # Write output matrix
    SeqIO.write(out_records, outfile, 'fasta')


if __name__ == '__main__':
     # Command line parsing via argparse
    parser = ap.ArgumentParser()
    parser.add_argument(type=str, dest='real_infile',
                        help='Path to fasta file with missing data to be replicated')
    parser.add_argument(type=str, dest='sim_infile',
                        help='Path to fasta file with missing data to be replicated')
    parser.add_argument(type=str, dest='out_file',
                        help='Path to output fasta file with replicated data')
    args = parser.parse_args()


    insert_missing_data(args.real_infile, args.sim_infile, args.out_file)
