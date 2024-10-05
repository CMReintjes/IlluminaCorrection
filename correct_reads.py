from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import argparse


def parse_arguments(test_args):
    "Parse command line arguments as input"
    parser = argparse.ArgumentParser(description='Error correction for Illumina sequencing data using k-mer counts')

    # Positional arguments
    parser.add_argument('file', type=str, help='Illumina reads file for correction')
    parser.add_argument('format', type=str, choices=['fasta', 'fastq'], help='Format of the Illumina file (.fasta or .fastq)')
    parser.add_argument('kmer_length', type=int, help='Length of the k-mer')
    parser.add_argument('threshold', type=int, help='Error threshold for correction')

    # Optional flag arguments
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')

    if test_args:
        return parser.parse_args(test_args)
    else:
        return parser.parse_args()


def read_illumina_file(args, file_name):
    "Parses the Illumina file and gets kmer frequency from sequences."
    file_format = args.format
    verbose = args.verbose
    k = args.kmer_length

    if verbose:
        print(f'Opening file: {file_name}')
    kmer_frequency = {}
    with open(file_name, 'r') as handle:
        if verbose:
            print(f'Parsing file as {file_format} format...')
        for record in SeqIO.parse(handle, file_format):
            if verbose:
                print(f'Processing record: {record.id}')
            get_kmer_frequency(str(record.seq), k, kmer_frequency)
    return kmer_frequency


def get_kmer_frequency(sequence, k, kmer_frequency):
    "Processes k-mers from a sequence and updates the frequency dictionary."
    kmer = sequence[:k]
    kmer_frequency[kmer] = kmer_frequency.get(kmer, 0) + 1
    for pos in range(1, len(sequence) - k + 1):
        kmer = kmer[1:] + sequence[pos + k - 1]
        kmer_frequency[kmer] = kmer_frequency.get(kmer, 0) + 1


def plot_kmer_frequencies(kmer_frequency, title):
    "Plots a histogram of k-mer frequencies."
    plt.figure(figsize=(10, 6))
    plt.hist(kmer_frequency.values(), bins=50)
    plt.xlabel('K-mer Frequency')
    plt.ylabel('Count')
    plt.title(title)
    plt.show()


def check_kmer_error(kmer, kmer_frequency, threshold):
    "Checks if a k-mer is below the error threshold. Returns True if it has an error."
    try:
        return kmer_frequency[kmer] < threshold
    except KeyError:
        return True  # If k-mer is not found, assume it has an error


def make_window_correction(args, window, kmer_frequency, threshold, forward):
    "Makes corrections to an erroneous base by trying all possible base changes."
    nucleotides = list('ATCG')
    if forward:
        original = window[-1]
        position = -1
    else:
        original = window[0]
        position = 0

    corrected_window = window[:]
    for n in nucleotides:
        if n != original:
            corrected_window[position] = n
            corrected_kmer = ''.join(corrected_window)

            # Verify if the corrected k-mer is above the error threshold
            if not check_kmer_error(corrected_kmer, kmer_frequency, threshold):
                if args.verbose:
                    print(f'Correction successful for k-mer: {corrected_kmer}')
                return corrected_window

    # If no correction was successful, revert the change
    if args.verbose:
        print(f'Correction failed for k-mer: {"".join(window)}')
    return window


def check_sequences(args, kmer_frequency):
    "Parses an Illumina sequencing file, and then corrects erroneous bases."
    file_path = args.file
    file_format = args.format
    verbose = args.verbose
    threshold = args.threshold
    k = args.kmer_length
    corrected_sequences = []

    if verbose:
        print(f'Opening file: {file_path}')
    with open(file_path, 'r') as input_file:
        if verbose:
            print(f'Parsing file as {file_format} format...')
        for record in SeqIO.parse(input_file, file_format):
            if verbose:
                print(f'Processing record: {record.id}')
            sequence = str(record.seq)
            seq_list = list(sequence)
            error_pos = -1

            # Forward pass through the sequence
            for pos in range(len(seq_list) - k + 1):
                kmer = ''.join(seq_list[pos:pos+k])
                if check_kmer_error(kmer, kmer_frequency, threshold):
                    if args.verbose:
                        print(f'K-mer {kmer} is below the error threshold.')
                    window = seq_list[pos:pos+k]
                    new_window = make_window_correction(args, 
                        window, kmer_frequency, threshold, forward=True)
                    seq_list[pos:pos+k] = list(new_window)
                    if args.verbose:
                        print(f'Corrected nucleotides: {new_window}')
                else:
                    if error_pos == -1:
                        error_pos = pos

            # Backward pass if the first k-mer has an error
            if error_pos > 0:
                for pos in range(error_pos - 1, -1, -1):
                    kmer = ''.join(seq_list[pos:pos+k])
                    if check_kmer_error(kmer, kmer_frequency, threshold):
                        if args.verbose:
                            print(f'Backward pass: K-mer {kmer} is below the error threshold.')
                        window = seq_list[pos:pos+k]
                        new_window = make_window_correction(args, 
                            window, kmer_frequency, threshold, forward=False)
                        seq_list[pos:pos+k] = list(new_window)
                        if args.verbose:
                            print(f'Backward pass corrected nucleotides: {new_window}')

            corrected_sequence = ''.join(seq_list)
            corrected_record = record
            corrected_record.seq = Seq(corrected_sequence)
            # Write corrected sequences to output file
            output_name = 'corrected_output.' + file_format
            with open(output_name, 'a') as output_file:
                SeqIO.write(corrected_record, output_file, file_format)


def main():
    test_args = ['reads_to_correct.fq', 'fastq', '21', '2', '-v']
    args = parse_arguments(test_args)
    if args.verbose:
        print(args)
    corrected_file_name = "corrected_output."+args.format
    with open(corrected_file_name, 'w'):
        pass
    kmer_frequency = read_illumina_file(args, args.file)
    plot_kmer_frequencies(kmer_frequency, 'K-mer Frequencies Before Correction')
    check_sequences(args, kmer_frequency)
    # Recalculate k-mer counts after correction
    corrected_kmer_frequency = read_illumina_file(args, corrected_file_name)
    plot_kmer_frequencies(corrected_kmer_frequency,
                          'K-mer Frequencies After Correction')


main()
