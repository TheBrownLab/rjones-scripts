#/bin/usr/env python3

import glob
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse as ap


def parse_input_proteome(input_proteome, unique_id):
    '''
    Parses input proteome file and returns a dictionary of the sequences and original headers.

    :param input_proteome: path to input proteome
    :type input_proteome: str
    :param unique_id: unique identifier for the proteome
    :type unique_id: str
    :return: dictionary of proteome contigs and their headers
    :rtype: dict
    '''
    data_dict = {}
    with open(input_proteome, 'r') as input_proteome:
        count = 0
        for record in SeqIO.parse(input_proteome, 'fasta'):
            data_dict[f'{unique_id}_g{count}'] = {
                'seq': str(record.seq),
                'header': record.description
            }
            count += 1
    
    return data_dict


def write_outputs(data_dict, output_proteome, key_tsv):
    '''
    Writes output proteome with new headers and key tsv files.

    :param data_dict: dictionary of proteome contigs and their headers
    :type data_dict: dict
    :param output_proteome: path to output proteome
    :type output_proteome: str
    :param key_tsv: path to key tsv file
    :type key_tsv: str
    '''
    # Initialize lists
    key_frags = []
    records = []

    # Append Seq records to records list
    for k, v in data_dict.items():
        key_frags.append(f'{k}\t{v["header"]}')
        records.append(SeqRecord(Seq(v['seq']),
                                 id=k,
                                 name='',
                                 description=''))

    # Write output proteome
    with open(output_proteome, 'w') as output_proteome:
        SeqIO.write(records, output_proteome, 'fasta')

    # Write key tsv
    with open(key_tsv, 'w') as output_tsv:
        output_tsv.write(f'new_header\told_header\n')
        output_tsv.write('\n'.join(key_frags))


if __name__ == '__main__':

    # Command line parsing via argparse
    parser = ap.ArgumentParser()
    parser.add_argument(type=str, dest='in_dir', help='Path to input dir containing proteomes.')
    parser.add_argument(type=str, dest='out_dir', help='Path to output directory.')
    args = parser.parse_args()

    # Make output dir
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)
    
    in_dir = os.path.abspath(args.in_dir)
    for file in glob.glob(f'{in_dir}/*.faa'):
        # Get full path to proteome and make sure it exists
        input_proteome = os.path.abspath(file)
        print(f'Processing {input_proteome}')

        # Use file basename as unique id
        unique_id = os.path.splitext(input_proteome)[0].split('/')[-1]

        # Get output proteome key names
        output_proteome = f'{out_dir}/{unique_id}.faa'
        key_tsv = f'{out_dir}/{unique_id}.tsv'

        # Parse input proteome
        data_dict = parse_input_proteome(input_proteome, unique_id)
        
        # Write output proteome and key
        write_outputs(data_dict, output_proteome, key_tsv)