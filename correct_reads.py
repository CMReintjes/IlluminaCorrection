#!/usr/bin/env python

"""
Identify kmers from Illumina Sequences data to correct errors within the dataset.
Plot the distribution of the kmers on a histogram before and after correction
"""

# General Imports
import argparse
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
       print(f'Creating dictionary of only trusted kmers using threshold of: {limit}')
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
                    kmer = record.seq[pos:pos+kmer_length]
                    kmer_frequency = get_Kmer_Frequency(
                        kmer, kmer_frequency)
    except FileNotFoundError:
        print(f'File {file} not found')
    if args.verbose:
       print('Successfully created dictionary of kmer frequencies')
    
    return kmer_frequency


def get_Hamming_Distance(kmer, trusted):
    "Calculate and return the Hamming Distance between the erroneous and trusted kmers"
    pass



def check_Sequence(args, kmer_trusted, sequence):
    "Check the number of kmers within the sequence that are below the threshold"
    errors = []
    kmer_length = args.kmer_size
    for pos in range(0, len(sequence)-kmer_length+1,int((kmer_length/3))):
        kmer = sequence[pos:pos+kmer_length]
        if kmer not in kmer_trusted.keys():
            errors.append(pos)

    if len(errors) == 0:
        return None
    return errors


def correct_Sequences(args, kmer_trusted):
    "Correct erroneous sequences within the file"
    file = args.file
    file_format = args.format
    kmer_length = args.kmer_size
    corrected_kmers = {}
    trusted_array = np.array([list(trusted) for trusted in kmer_trusted.keys()])
    # print(trusted_array)
    '''with open('errors.csv', 'w') as output_file:
        pass'''
    with open('corrected.fastq', mode='w') as correctFile:
        pass

    # Iterate through sequence file
    with open(file, mode='rt') as seqFile:
        sio = SeqIO.parse(seqFile, file_format)
        count = 0
        for record in sio:
            count += 1
            sequence = str(record.seq)
            errors = check_Sequence(args, kmer_trusted, sequence)
            if errors is not None:
                '''with open('errors.csv', 'a') as output_file:
                    errors.sort()
                    output_file.write(f'\n{sequence}')
                    [output_file.write(f',{errors[i]}') for i in range(len(errors))]'''
                error_kmers = [sequence[pos:pos+kmer_length] for pos in errors]
                print(error_kmers)

                for error_kmer in error_kmers:
                    errors = np.array(error_kmer)
                    distances = np.sum(errors != trusted_array, axis=1)
                    index = np.argmin(distances)
                    closest_kmer = ''.join(list(trusted_array)[index])
                    if error_kmer not in corrected_kmers:
                        corrected_kmers[error_kmer] = closest_kmer
                
                corrected_sequence = sequence
                for pos in range(len(corrected_sequence)):
                    kmer = corrected_sequence[pos:pos+kmer_length]
                    if kmer in kmer_trusted and kmer in corrected_kmers:
                        corrected_sequence = corrected_sequence[:pos] + corrected_kmers[kmer] + corrected_sequence[pos+kmer_length:]
                print(f'{count}: {corrected_sequence}')
                with open('corrected.fastq', mode='a') as correctFile:
                    record.seq = Seq(corrected_sequence)
                    SeqIO.write(record, correctFile, file_format)
                        #print_Sequences(sequence, pos, kmer, corrected_kmers[kmer])
                    


def main(args):
    "Main Function"
    corrected = 'corrected.'+args.format
    kmer_frequency = get_Sequences(args, file=args.file)
    kmer_trusted = get_Trusted(args, kmer_frequency)
    plot_Histogram(kmer_frequency)
    correct_Sequences(args, kmer_trusted)
    plot_Histogram(get_Sequences(args, corrected))


args = parse_Args()
main(args)
