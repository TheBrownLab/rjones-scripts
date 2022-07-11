#!/bin/usr/env python3

import os
import argparse as ap

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC


def plot_gc(input_genome, output_plot):
    '''
    Plots GC content of input genome.

    :param input_genome: path to input genome
    :type input_genome: str
    '''
    gc_content = []
    with open(input_genome, 'r') as input_genome:
        for record in SeqIO.parse(input_genome, 'fasta'):
                gc_content.append(GC(record.seq))
    
    plt.hist(gc_content, bins=100)
    plt.savefig(output_plot)
    plt.close()


def filter_by_gc(input_genome, output_genome, gc_threshold):
    '''
    Filters input genome by GC content.

    :param input_genome: path to input genome
    :type input_genome: str
    :param output_genome: path to output genome
    :type output_genome: str
    :param gc_threshold: GC content threshold
    :type gc_threshold: float
    '''
    with open(input_genome, 'r') as input_genome, open(output_genome, 'w') as output_genome:
        for record in SeqIO.parse(input_genome, 'fasta'):
            if GC(record.seq) < gc_threshold:
                SeqIO.write(record, output_genome, 'fasta')



if __name__ == '__main__':
     # Command line parsing via argparse
    parser = ap.ArgumentParser()
    parser.add_argument(type=str, dest='in_file', help='Path to input genome')
    parser.add_argument(type=str, dest='out_dir', help='Path to output directory')
    parser.add_argument('--filter', type=float, dest='filter', help='Output filtered genome with GC content below this value.', default=None)
    parser.add_argument('--plot', dest='plot', action='store_true', help='Plot GC content histogram')
    args = parser.parse_args()
    
    # Make output dir
    out_dir = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    # Input genome path
    in_file = os.path.abspath(args.in_file)

    # Output genome path
    basename = os.path.splitext(args.in_file)[0].split('/')[-1]
    output_genome = os.path.join(args.out_dir, f'{basename}.gc_{args.filter}_filtered.fas')
    output_plot = os.path.join(args.out_dir, f'{basename}.gc_histogram.png')


    if args.plot:
        plot_gc(in_file, out_dir)
    if args.filter:
        filter_by_gc(in_file, output_genome, args.filter)
