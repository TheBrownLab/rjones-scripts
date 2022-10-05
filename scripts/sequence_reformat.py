'''
Reformats a sequence file to a different format.
'''

import argparse as ap
from Bio import SeqIO

# Accepted out formats with respective suffix
EXT_DICT = {'fasta':          'fas',
            'phylip':         'phy',
            'phylip-relaxed': 'phy',
            'nexus':          'nex'}

if __name__ == '__main__':
    # Command line parsing via argparse
    parser = ap.ArgumentParser()
    parser.add_argument(type=str, dest='in_file', help='Path to input sequence file')
    parser.add_argument(type=str, dest='in_format', help='Input sequence format')
    parser.add_argument(type=str, dest='out_file', help='Path to output sequence file')
    parser.add_argument(type=str, dest='out_format', help='Output sequence format')
    args = parser.parse_args()

    # Initialize records list
    records = []

    # Parse input sequence file
    with open(args.in_file, 'r') as infile:
        for record in SeqIO.parse(infile, args.in_format):
            records.append(record)

    # Write output sequence file
    with open(args.out_file, 'w') as outfile:
        SeqIO.write(records, outfile, args.out_format)
