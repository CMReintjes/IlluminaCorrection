#!/usr/bin/env python

"""
Identify kmers from Illumina Sequences data to correct errors within the dataset.
Plot the distribution of the kmers on a histogram before and after correction
"""

# General Imports
import argparse
from cmath import e
from os import error
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Biopython Imports
import Bio
from Bio import SeqIO
from Bio.Seq import Seq


def parse_Args():
    "Initialize argument parser"
    "correct_reads.py <reads to correct> <format> <kmer size> <threshold>"
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


def plot_Histogram(frequency):
    "Plot the kmer frequency using a seaborn histogram"
    if args.verbose:
        print('Creating Histogram of kmer frequency')

    sns.histplot(list(frequency.values()), bins=50)
    plt.xlabel('Kmer Frequency')
    plt.ylabel('Number of Kmers')
    plt.title('Histogram of Kmer Frequency Distribution Values')
    plt.show()


def print_Sequences(sequence, pos, kmer, new_kmer):
    "Print out the new sequences and replacements when verbose is active"
    if args.verbose:
        print('Correcting sequence')
        [print(' ', end='') for j in range(pos)]
        print(kmer)
        [print(' ', end='') for j in range(pos)]
        print(new_kmer)
        print(sequence)


def append_Output(args, record):
    # Append the corrected sequences and sequence objects to the output file
    outputName = args.output+'.'+args.format
    if args.verbose:
        # print(f'Attempting to open output file: {outputName}')
        pass
    try:
        with open(file=outputName, mode='a') as seqOut:
            SeqIO.write(record, seqOut, args.format)
            if args.verbose:
                print(f'Successfully wrote record to {outputName}')
    except FileNotFoundError:
        if args.verbose:
            print('Output file not found. Creating file and opening')
        with open(file=outputName, mode='w') as seqOut:
            if args.verbose:
                print(f'Output file creation successful')
            SeqIO.write(record, seqOut, args.format)
            if args.verbose:
                print(f'Successfully wrote record to {outputName}')


def get_Kmer_Frequency(kmer, kmer_frequency):
    # Iterate through the file and hash kmers
    # If kmer in hash, increment frequency, if not, add key with value 1
    if kmer in kmer_frequency:
        kmer_frequency[kmer] += 1
    else:
        kmer_frequency[kmer] = 1

    return kmer_frequency


def get_Trusted(args, all_frequencies):
    "Return a dictionary containing only the trusted sequences based on a cutoff"
    limit = args.threshold + 1
    if args.verbose:
        print(
            f'Creating dictionary of only trusted kmers using threshold of: {limit}')
    trusted = {}
    for kmer, freq in all_frequencies.items():
        if freq > limit:
            trusted[kmer] = freq
    if args.verbose:
        print('Successfully created dictionary of trusted kmers')

    return trusted


def get_Sequences(args, file):
    # Get the frequency at which each kmer appears within the illumina file
    kmer_frequency = {}
    file_format = args.format
    kmer_length = args.kmer_size

    if args.verbose:
        print(f'Calculating frequency of all kmers within: {file}')

    try:
        # Iterate through sequence file
        with open(file, mode='rt') as seqFile:
            if args.verbose:
                print(f'Successfully opened file: {file}')
            sio = SeqIO.parse(seqFile, file_format)
            if args.verbose:
                print(f'Successfully parsed file with format: {file_format}')
            # Iterate through position in nucleotide sequence
            for record in sio:
                # print(record.seq)
                for pos in range(len(record.seq)-kmer_length+1):
                    # Create sub sequence of length n using string position
                    kmer = str(record.seq[pos:pos+kmer_length])
                    kmer_frequency = get_Kmer_Frequency(
                        kmer, kmer_frequency)
    except FileNotFoundError:
        print(f'File {file} not found')
    if args.verbose:
        print('Successfully created dictionary of kmer frequencies')

    return kmer_frequency


def correct_Kmers(args, pos, sequence, trusted_kmers):
    kmer_length = args.kmer_size
    nucleotides = ['A', 'C', 'G', 'T']
    trusted = set(trusted_kmers.keys())
    "Correct and return the erroneous kmers within the sequence"
    checking_kmers = []
    for i in range(kmer_length):
        check_list = list(sequence[pos:pos+kmer_length])
        options = [check_list for j in range(4)]
        for k in range(4):
            options[k][(kmer_length-i)-1] = nucleotides[k]
            #print(f'Option: {''.join(options[k])}')
            #print(f'Sequence: {options[k]}, Nucleotide: {nucleotides[k]}')

            if ''.join(options[k]) in trusted_kmers:
                print("kmers exist within the dictionary")

        




def check_Sequence(args, kmer_trusted, sequence):
    "Check the number of kmers within the sequence that are below the threshold"
    errors = []
    reverse = []
    initial_error = False
    initial_pos = 0
    needed_reverse = False
    kmer_length = args.kmer_size
    for pos in range(len(sequence)-kmer_length+1):
        kmer = sequence[pos:pos+kmer_length]
        if kmer not in kmer_trusted.keys() and pos == 0:
            initial_error = True
            needed_reverse = True
        elif kmer not in kmer_trusted.keys() and initial_error is True:
            continue
        elif kmer in kmer_trusted.keys() and initial_error is True:
            initial_error = False
            initial_pos = pos
            continue
        elif kmer not in kmer_trusted.keys() and initial_error is False:
            errors.append(pos)
        
    if needed_reverse:
        for pos in range(initial_pos, -1, -1):
            if kmer not in kmer_trusted.keys():
                reverse.append(pos)
    
    if len(errors) == 0 and len(reverse) == 0:
        return None, None
    elif len(errors) == 0:
        return None, reverse
    elif len(reverse) == 0:
        return errors, None
    else:
        return errors, reverse


def correct_Sequences(args, kmer_trusted):
    "Correct erroneous sequences within the file"
    file = args.file
    file_format = args.format
    kmer_length = args.kmer_size
    corrected_kmers = {}
    # print(trusted_array)
    with open('errors.csv', 'w') as output_file:
        pass
    with open('corrected.fastq', mode='w') as correctFile:
        pass

    # Iterate through sequence file
    with open(file, mode='rt') as seqFile:
        sio = SeqIO.parse(seqFile, file_format)
        count = 0
        for record in sio:
            count += 1
            sequence = str(record.seq)
            errors, reverse = check_Sequence(args, kmer_trusted, sequence)
            if errors is not None or reverse is not None:
                with open('errors.csv', 'a') as output_file:
                    output_file.write(f'\n{sequence}')
                    if errors is not None:
                        errors.sort()
                        [output_file.write(f',{errors[i]}') for i in range(len(errors))]
                        last = 0
                        # error_kmers = [sequence[pos:pos+kmer_length] for pos in errors]
                        for pos in errors:
                            correct_Kmers(args, pos, sequence, kmer_trusted)
                            break
                            
                    if reverse is not None:
                        reverse.sort()
                        [output_file.write(f',{reverse[i]}') for i in range(len(reverse))]
                break
                



def main(args):
    "Main Function"
    corrected = 'corrected.'+args.format
    kmer_frequency = get_Sequences(args, file=args.file)
    kmer_trusted = get_Trusted(args, kmer_frequency)
    #plot_Histogram(kmer_frequency)
    correct_Sequences(args, kmer_trusted)
    # plot_Histogram(get_Sequences(args, corrected))


args = parse_Args()
main(args)
