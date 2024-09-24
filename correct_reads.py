#!/usr/bin/env python

"""
Identify kmers from Illumina Sequences data to correct errors within the dataset.
Plot the distribution of the kmers on a histogram before and after correction
"""

# General Imports
import argparse
from typing import final
import numpy.random as npr

# Biopython Imports
from Bio import SeqIO
from Bio import Seq


def parseArgs():
    # Initialize argument parser
    # correct_reads.py <reads to correct> <format> <kmer size> <threshold>
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
    try:
        # Iterate through sequence file
        with open(file=args.file, mode='rt') as seqFile:
            if args.verbose:
                print(f'Successfully opened file: {args.file}')
            sio = SeqIO.parse(seqFile, args.format)
            if args.verbose:
                print(f'Successfully parsed file with format: {args.format}')
            # Iterate through position in nucleotide sequence
            for record in sio:
                # print(record.seq)
                for pos in range(len(record.seq)-kmer_length+1):
                    # Create sub sequence of length n using string position
                    kmer = record.seq[pos:pos+kmer_length]
                    # [print(' ',end='') for i in range(pos)]
                    # print(kmer)
                    kmer_frequency = getKmerFrequency(kmer, kmer_frequency)
    except FileNotFoundError:
        print(f'File {args.file} not found')

    '''with open('sequence_output.txt', 'w') as outputFile:
        outputFile.write(str(kmer_frequency))'''
    return kmer_frequency


def changeNucleotide(option, reverse, nuc):
    # Change the nucleotide at position before returning sequence
    kmer_length = args.kmer_size
    nucleotides = ['A', 'C', 'G', 'T']
    if reverse != True:
        option[kmer_length -
                1] = nucleotides[nuc]
    else:
        option[0] = nucleotides[nuc]
    # Join kmer together
    option = ''.join(option)

    # Compare all variants in <options> to determine error free kmer
    return option


def checkKmerCounts(kmer, kmer_frequency, reverse):
    # Check the kmer against the count within the hashmap
    max_kmer = None
    # Make list of kmer options
    options = [list(kmer) for n in range(4)]
    # print(f'Options:')
    for i in range(4):
        options[i] = changeNucleotide(option=options[i], reverse=reverse, nuc=i)
        if args.verbose:
            try:
                print(f'{options[i]}: {kmer_frequency[options[i]]}')
            except KeyError:
                print(f'{options[i]}: 0')
        # Compare all variants in <options> to determine error free kmer
        if options[i] in kmer_frequency:
            if max_kmer == None or kmer_frequency[options[i]] > kmer_frequency[max_kmer]:
                # print(f'{kmer_frequency[options[i]]}')
                max_kmer = options[i]
            # Replace kmer with the kmer with highest frequency
            kmer = max_kmer
    # Return kmer
    return kmer


def printSequences(kmer_length, sequence, pos, kmer, new_kmer):
    # Print out the new sequences and replacements when verbose is active
    pos-1
    if args.verbose:
        print('Correcting sequence')
        [print(' ', end='') for j in range(pos)]
        print(kmer)
        [print(' ', end='') for j in range(pos)]
        print(new_kmer)
        print(sequence)


def checkErrorState(frequency, threshold, pos, initialError):
    # Check the error state of the kmer based on threshold
    # Ensure that if there was an initial error, ignore until good kmer
    if frequency <= threshold and pos == 0:
        initialError == True
        return False, initialError
    elif frequency <= threshold and initialError == True:
        return False, initialError
    elif frequency > threshold and initialError == True:
        initialError = False
        return False, initialError
    elif frequency <= threshold:
        return True, initialError
    else:
        return False, False


def makeNewSequence(kmer_length, sequence, pos, kmer, newKmer):
    if kmer != newKmer:
        newSeq = str(sequence[:pos])+newKmer + \
            str(sequence[pos+kmer_length:])
        sequence = newSeq
        printSequences(kmer_length, sequence, pos, kmer, newKmer)
    return sequence


def appendOutput(record):
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


def checkSequences(kmer_frequency):
    # Check the sequences for errors using kmer counts
    kmer_length = args.kmer_size
    threshold = args.threshold
    try:
        with open(file=args.file, mode='rt') as seqFile:
            for record in SeqIO.parse(seqFile, args.format):
                sequence = record.seq
                initialError, error = False, False
                # Start reads at the beginning of sequence
                for pos in range(0, len(sequence)-kmer_length+1):
                    kmer = str(sequence[pos:pos+kmer_length])
                    # kmer = str(record.seq[pos:pos+kmer_length])
                    # Check if the kmer is under the error threshold and if there was a previous error
                    try:
                        error, initialError = checkErrorState(
                            kmer_frequency[kmer], threshold, pos, initialError)
                    except KeyError:
                        if args.verbose:
                            print("Forward Checking Sequence\nProblem fixing error: Kmer not in hash\n")
                            printSequences(kmer_length, sequence,
                                           pos, kmer, str(record.seq[pos:pos+kmer_length]))
                    if error:
                        # Perform check on kmer variants and return valid
                        newKmer = checkKmerCounts(
                            kmer, kmer_frequency, reverse=False)
                        sequence = makeNewSequence(
                            kmer_length, sequence, pos, kmer, newKmer)
                # Run error checking in reverse do fix initial errors and verify
                if args.verbose:
                    print('Reverse Checking Sequence')

                for pos in range(len(sequence)-kmer_length, -1, -1):
                    kmer = str(sequence[pos:pos+kmer_length])
                    # kmer = str(record.seq[pos:pos+kmer_length])
                    # [print(f'{key}:{value}') for key, value in kmer_frequency.items()]
                    if kmer_frequency[kmer] <= threshold:
                        newKmer = checkKmerCounts(
                            kmer, kmer_frequency, reverse=True)
                        sequence = makeNewSequence(
                            kmer_length, sequence, pos, kmer, newKmer)
                appendOutput(record)
    except FileNotFoundError:
        print(f'File {args.file} not found.')


def main():
    kmer_frequency = getSequences()
    # [print(f'{key}:{value}') for key, value in kmer_frequency.items()]
    checkSequences(kmer_frequency)


args = parseArgs()
main()
