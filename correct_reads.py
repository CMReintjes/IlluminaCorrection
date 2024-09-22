#!/usr/bin/env python

"""
Identify kmers from Illumina Sequences data to correct errors within the dataset.
Plot the distribution of the kmers on a histogram before and after correction
"""

# Imports
import Bio
import argparse
import sys
import numpy.random as npr
from Bio import SeqIO
from Bio import Seq

# Initialize argument parser

# correct_reads.py <reads to correct> <format> <kmer size> <threshold>


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str,
                        help='A file of Illumina reads that are to be corrected')
    parser.add_argument('format', type=str, choices=[
                        'fasta', 'fastq'], help='The sequence format of the reads to correct. Only FASTQ and FASTA are valid entries')
    parser.add_argument('kmer_size', type=int,
                        help='The kmer size to use for the analysis')
    parser.add_argument('threshold', type=int,
                        help='The minimum count a kmer must have to not be considered an error')

    args = parser.parse_args()
    return args


def getKmerFrequency(kmer, kmer_frequency):
    # Iterate through the file and hash kmers
    # If kmer in hash, increment frequency, if not, add key with value 1
    if kmer in kmer_frequency:
        kmer_frequency[kmer] += 1
    else:
        kmer_frequency[kmer] = 1

    return kmer_frequency


def getSequences():
    # Get the frequency at which each kmer appears within the illumina file
    kmer_frequency = {}
    kmer_length = args.kmer_size
    # file_format = args.format
    # Iterate through sequence file
    with open(file=args.file, mode='rt') as seqFile:
        sio = SeqIO.parse(seqFile, args.format)
        # Iterate through position in nucleotide sequence
        for record in sio:
            for pos in range(len(record.seq)-kmer_length):
                # Create sub sequence of length n using string position
                kmer = record.seq[pos:pos+kmer_length]
                kmer_frequency = getKmerFrequency(kmer, kmer_frequency)

    '''with open('sequence_output.txt', 'w') as outputFile:
        outputFile.write(str(kmer_frequency))'''
    return kmer_frequency


def changeNucleotide(option, reverse):
    # Change the nucleotide at position before returning sequence
    kmer_length = args.kmer_size
    if reverse != True:
        option[kmer_length -
               1] = npr.choice(list('ATGC'), replace=False)
    else:
        option[0] = npr.choice(list('ATGC'), replace=False)
    # Join kmer together
    option = ''.join(option)
    # Compare all variants in <options> to determine error free kmer
    return option


def checkKmerCounts(kmer, kmer_frequency, reverse):
    # Check the kmer against the count within the hashmap
    threshold = args.threshold
    max_kmer = None
    # Check if the kmer is under the error threshold
    if kmer_frequency[kmer] <= threshold:
        # Make list of kmer options
        options = [list(kmer) for n in range(4)]
        # print(f'Options:')
        for i in range(4):
            options[i] = changeNucleotide(option=options[i], reverse=reverse)
            # Compare all variants in <options> to determine error free kmer
            if options[i] in kmer_frequency:
                # print(f'{options[i]}: ', end='')
                if max_kmer == None or kmer_frequency[options[i]] > kmer_frequency[max_kmer]:
                    # print(f'{kmer_frequency[options[i]]}')
                    max_kmer = options[i]
                # Replace kmer with the kmer with highest frequency
                kmer = max_kmer
    # Return kmer
    return kmer


def checkSequences(kmer_frequency):
    # Check the sequences for errors using kmer counts
    kmer_length = args.kmer_size
    with open(file=args.file, mode='rt') as seqFile:
        for record in SeqIO.parse(seqFile, args.format):
            sequence = record.seq
            errors = True
            for pos in range(kmer_length, len(sequence)-kmer_length):
                kmer = record.seq[pos:pos+kmer_length]
                newKmer = checkKmerCounts(kmer, kmer_frequency, reverse=False)
                if kmer != newKmer:
                    newSeq = str(sequence[:pos])+newKmer+str(sequence[pos+kmer_length:])
                    sequence = newSeq
            for pos in range(kmer_length, -1, -1):
                kmer = record.seq[pos:pos+kmer_length]
                newKmer = checkKmerCounts(kmer, kmer_frequency, reverse=True)
                if kmer != newKmer:
                    newSeq = str(sequence[:pos])+newKmer+str(sequence[pos+kmer_length:])
                    sequence = newSeq


def main():
    kmer_frequency = getSequences()
    checkSequences(kmer_frequency)


args = parseArgs()
main()
