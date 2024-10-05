#!/usr/bin/env python

"""
Identify kmers from Illumina Sequences data to correct errors within the dataset.
Plot the distribution of the kmers on a histogram before and after correction
"""

# General Imports
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Biopython Imports
from Bio import SeqIO

def parse_args():
    "Initialize argument parser correct_reads.py <reads to correct> <format> <kmer size> <threshold>"
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str,
                        help='A file of Illumina reads that are to be corrected')
    parser.add_argument('format', type=str, choices=[
                        'fasta', 'fastq'], help='The sequence format of the reads to correct. Only FASTQ and FASTA are valid entries')
    parser.add_argument('kmer_size', type=int,
                        help='The kmer size to use for the analysis')
    parser.add_argument('threshold', type=int,
                        help='The minimum count a kmer must have to not be considered an error')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('-o', '--output', type=str, default='corrected_reads',
                        help='Enable secondary output file name')

    args = parser.parse_args()
    return args


def plot_histogram(frequency, bins):
    "Plot the kmer frequency using a seaborn histogram"
    sns.histplot(list(frequency.values()), bins=bins)
    plt.xlabel('Kmer Frequency')
    plt.ylabel('Number of Kmers')
    plt.title('Histogram of Kmer Frequency Distribution Values')
    plt.show()


def print_sequences(args, sequence, pos, kmer, new_kmer):
    "Print out the new sequences and replacements when verbose is active"
    if args.verbose:
        print('Correcting sequence')
        spacing = ''.join([' ' for j in range(pos)])
        print(f'{spacing}{kmer}')
        print(f'{spacing}{new_kmer}')
        print(sequence)


def append_output(args, record):
    "Append the corrected sequences and sequence objects to the output file"
    output_name = args.output+'.'+args.format
    if args.verbose:
        # print(f'Attempting to open output file: {outputName}')
        pass
    try:
        with open(output_name,'a') as seq_out:
            SeqIO.write(record, seq_out, args.format)
            if args.verbose:
                pass
                #print(f'Successfully wrote record to {output_name}')
    except FileNotFoundError:
        if args.verbose:
            print('Output file not found. Creating file and opening')
        with open(output_name,'w') as seq_out:
            if args.verbose:
                print('Output file creation successful')
            SeqIO.write(record, seq_out, args.format)
            if args.verbose:
                print(f'Successfully wrote record to {output_name}')


def check_error_state(frequency, threshold, pos, initial_error):
    "Check the error state of the kmer based on threshold. Ensure that if there was an initial error, ignore until good kmer"
    if frequency <= threshold and not pos:
        initial_error = True
        return False, initial_error
    elif frequency <= threshold and initial_error is True:
        return False, initial_error
    elif frequency > threshold and initial_error is True:
        initial_error = False
        return False, initial_error
    elif frequency <= threshold:
        return True, initial_error
    else:
        return False, False


def get_kmer_frequency(kmer, kmer_frequency):
    "Iterate through the file and hash kmers. If kmer in hash, increment frequency, if not, add key with value 1"
    if kmer in kmer_frequency:
        kmer_frequency[kmer] += 1
    else:
        kmer_frequency[kmer] = 1

    return kmer_frequency


def get_sequences(args, file):
    "Get the frequency at which each kmer appears within the illumina file"
    kmer_frequency = {}
    file_format = args.format
    kmer_length = args.kmer_size
    try:
        # Iterate through sequence file
        with open(file, 'rt') as seq_file:
            if args.verbose:
                print(f'Successfully opened file: {file}')
            sio = SeqIO.parse(seq_file, file_format)
            if args.verbose:
                print(f'Successfully parsed file with format: {file_format}')
            # Iterate through position in nucleotide sequence
            for record in sio:
                # print(record.seq)
                for pos in range(len(record.seq)-kmer_length+1):
                    # Create sub sequence of length n using string position
                    kmer = record.seq[pos:pos+kmer_length]
                    kmer_frequency = get_kmer_frequency(
                        kmer, kmer_frequency)
    except FileNotFoundError:
        print(f'File {file} not found')
    return kmer_frequency


def check_sequences(args, kmer_frequency):
    "Check the frequency of kmers against the hash for errors, and correct erroneous bases"
    file_name = args.file
    file_format = args.format
    kmer_length = args.kmer_size
    threshold = args.threshold
    with open(file_name, 'rt') as seq_file:
        if args.verbose:
            print(f'Successfully opened file: {file_name}')
        sio = SeqIO.parse(seq_file, file_format)
        for r in sio:
            sequence = list(r.seq)
            for n in range(len(sequence)-kmer_length+1):
                kmer = ''.join(sequence[n:n+kmer_length])
                if kmer_frequency[kmer] < threshold:
                    nucleotides = list('ATGC')
                    for i in nucleotides:
                        option = list(kmer)
                        option[kmer_length-1] = i
                        new_kmer = ''.join(option)
                        if new_kmer in kmer_frequency and kmer_frequency[new_kmer] >= threshold:
                            sequence[n+kmer_length] = i
                            break
            print(''.join(sequence))
                        

def main():
    "Main Function"
    args = parse_args()
    kmer_frequency = get_sequences(args, file=args.file)
    #plot_histogram(kmer_frequency, 50)
    check_sequences(args, kmer_frequency)
    adjusted = get_sequences(args, 'corrected_reads.'+args.format)
    #plot_histogram(adjusted, 50)


main()
