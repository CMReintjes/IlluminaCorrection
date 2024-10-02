#!/usr/bin/env python

"""
Identify kmers from Illumina Sequences data to correct errors within the dataset.
Plot the distribution of the kmers on a histogram before and after correction
Current command: py correct_reads.py reads_to_correct.fq fastq 21 2 -v
"""

# General Imports
import argparse
import seaborn as sns
import matplotlib.pyplot as plt

# Biopython Imports
from Bio import SeqIO
from Bio.Seq import Seq

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
    "Print out the new sequences and replacements when verbose is active (OLD CODE, CAN REMOVE)"
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
    try: # Error handling for file not found
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
    if frequency <= threshold and not pos: # Check if first kmer has error
        initial_error = True
        return False, initial_error
    elif frequency <= threshold and initial_error is True: # Ignore following kmers after initial error
        return False, initial_error
    elif frequency > threshold and initial_error is True: # Check for good kmer following initial error
        initial_error = False
        return False, initial_error
    elif frequency <= threshold: # Once good is found, start correcting erroneous kmers
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
    try: # Error handling for missing file
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
                    kmer = str(record.seq[pos:pos+kmer_length])
                    kmer_frequency = get_kmer_frequency(
                        kmer, kmer_frequency)
    except FileNotFoundError:
        print(f'File {file} not found')
    return kmer_frequency


def correct_sequence(args, kmer_slots):
    "Gets corrected sequence from corrected kmer slots"
    sequence = ''
    for kmer in kmer_slots:
        sequence += kmer[0] # Concatenante strings with the first base of every kmer
    sequence += kmer[1:] # Add the last bases from final kmer at end

    if args.verbose:
        print(print(f'Corrected Sequence:\n{sequence}'))

    return sequence


def make_new_sequence(args, sequence, pos, kmer, new_kmer):
    "Make new sequence using the adjusted kmer and sequence (OLD CODE, NOT NEEDED)"
    kmer_length = args.kmer_size
    new_seq = str(sequence[:pos])+new_kmer + \
        str(sequence[pos+kmer_length:])
    sequence = new_seq
    print_sequences(args, sequence, pos, kmer, new_kmer)
    return sequence


def check_kmer_counts(args, pos, kmer, kmer_slots, kmer_frequency, reverse):
    "Check the kmer against the count within the hashmap"
    kmer_length = args.kmer_size
    max_kmer = kmer
    # Make list of kmer options
    options = [list(kmer) for n in range(4)]
    nucleotides = ['A', 'C', 'G', 'T']
    max_test = list(kmer_slots) # re-declare list since immulatable
    # Change nucleotide at position within options
    for i in range(4):
        kmer_test = list(kmer_slots)
        # Start the for loop at current kmer position
        future_start = pos
        if reverse is True:
            options[i][0] = nucleotides[i] # Change nucleotide at first position
            point = 0 # Point = position of nucleotide change within kmer
            future_stop = pos-kmer_length-1 # Future check stop point
            future_step = -1 # Future point step, -1 since backwards

        else:
            options[i][kmer_length-1] = nucleotides[i]
            point = kmer_length -1
            future_stop = pos+kmer_length+1
            future_step = 1
        # Check if this variant is in the hash
        option = ''.join(options[i])
        if option in kmer_frequency:
            future = False
            # If kmer is in hash, iterate through future kmers and change nucleotide to match
            for j in range(future_start, future_stop, future_step):
                if j >= len(kmer_test) or j >= kmer_length: # prevent index errors by breaking out of loop at kmer end points
                    break
                test = list(kmer_slots[j]) # convert string to list of bases
                test[point] = nucleotides[i] # change base at point
                kmer_test[j] = ''.join(test) # Join list back to string for hash check
                # Update point for next kmer
                if reverse:
                    point +=1
                else:
                    point -=1

                # Check if future kmer in hash
                if kmer_test[j] not in kmer_frequency:
                    future = False
                    break

                # If reaching end, future kmers exist
                if j == future_stop or j == len(kmer_test):
                    future = True

            # Update max with current options that exist
            if future is True and kmer_frequency[option] >= kmer_frequency[kmer]:
                # print(f'{kmer_frequency[options[i]]}')
                max_kmer = option
                max_test = list(kmer_test)
                future = False
                #print(f'Updated kmer_slots with: {kmer_test}')
    # Return kmer and kmer list
    return max_kmer, max_test


def check_sequences(args, kmer_frequency):
    "Check the sequences for errors using the kmer counts"
    file_name = args.file
    file_format = args.format
    kmer_length = args.kmer_size
    threshold = args.threshold

    with open(file_name, 'rt') as seq_file:
        for record in SeqIO.parse(seq_file, file_format):
            sequence = str(record.seq) # Create a sequence from seq object
            error, initial_error = False, False # Bools for initial error handling
            start_error = 0 # Tracker for initial error position
            if args.verbose: # verbose output flag
                print(f'{record.id}:\n{sequence}')

            # Starts reads at the beginning of the sequence and put in a list
            kmer_slots = [None for i in range(len(sequence)-kmer_length+1)]
            for pos in range(0, len(kmer_slots)):
                kmer_slots[pos] = str(sequence[pos:pos+kmer_length])

            for pos in range(len(kmer_slots)):
                kmer = kmer_slots[pos]
                # print(f'Checking kmer: {kmer} at position {pos}')
                # Check if the kmer is under the error threshold and if there was a previous error
                error, initial_error = check_error_state(
                kmer_frequency[kmer], threshold, pos, initial_error)
                if initial_error:
                    start_error = pos+kmer_length+1 # Update initial error position tracker
                if error:
                    new_kmer, kmer_slots = check_kmer_counts(args, pos, kmer, kmer_slots, kmer_frequency, reverse=False)
                    record.seq = Seq(correct_sequence(args, kmer_slots)) # Update sequence in record

            # If there was an initial error recheck start of sequence
            if start_error > 0:
                if args.verbose:
                    print(f'Reverse reading sequence at {start_error} due to initial error.')
                # Start at known good kmer position and work backwards
                for rev in range(start_error, -1, -1):
                    kmer = kmer_slots[rev]
                    if kmer_frequency[kmer] <= threshold:
                        new_kmer, kmer_slots = check_kmer_counts(
                            args, pos, kmer, kmer_slots, kmer_frequency, reverse=True)
                        record.seq = Seq(correct_sequence(args, kmer_slots))
            append_output(args, record)


def main():
    "Main Function"
    args = parse_args()
    correct_file_name = 'corrected_reads.'+args.format
    with open(correct_file_name, 'w'):
        pass

    kmer_frequency = get_sequences(args, file=args.file)
    plot_histogram(kmer_frequency, 50)
    check_sequences(args, kmer_frequency)
    adjusted = get_sequences(args, correct_file_name)
    plot_histogram(adjusted, 50)


main()
