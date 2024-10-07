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
    kmer_frequency = {}

    if verbose:
        print(f'Opening file: {file_name}')
    # Open Illumina File
    with open(file_name, 'r') as input_file:
        if verbose: # Verbose output
            print(f'Parsing file as {file_format} format...')
        # Iterate through records
        for record in SeqIO.parse(input_file, file_format):
            if verbose: # Verbose output
                print(f'Processing record: {record.id}')
            # Get kmer frequency using dictionary
            kmer_frequency = get_kmer_frequency(str(record.seq), k, kmer_frequency)
    return kmer_frequency


def get_kmer_frequency(sequence, k, kmer_frequency):
    "Gets k-mers from the sequence and updates the frequency dictionary."
    for pos in range(len(sequence)-k+1):
        # Create kmers of length k using string positions
        kmer = sequence[pos:pos+k]
        if kmer in kmer_frequency: # If key exists add one to frequency
            kmer_frequency[kmer] += 1
        else: # If new kmer then add with frequency 1
            kmer_frequency[kmer] = 1

    return kmer_frequency


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
    try: # Check if kmer is below the error threshold
        if kmer_frequency[kmer] < threshold:
            return True
    except KeyError: # Error handling just in case
        return True  


def make_window_correction(args, window, kmer_frequency, threshold, forward):
    "Makes corrections to an erroneous base by trying all possible base changes."
    nucleotides = list('ATCG') # List of nucleotides for base changes
    if forward: # If checking in forward direction use last base
        position = -1
        original = window[position]
        
    else: # If checking in reverse direction use first base
        position = 0
        original = window[position]
        
    corrected_window = list(window) # Make copy of original
    for nuc in nucleotides:
        if nuc != original:
            corrected_window[position] = nuc
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

    if verbose:
        print(f'Opening file: {file_path}')
    # Open file
    with open(file_path, 'r') as input_file:
        if verbose:
            print(f'Parsing file as {file_format} format...')
        # Iterate through the records
        for record in SeqIO.parse(input_file, file_format):
            if verbose:
                print(f'Processing record: {record.id}')
            # Create sequences based on the record objects
            sequence = str(record.seq)
            seq_list = list(sequence) # Convert to list
            error_pos = -1 # For keeping track of error position if there is initial error

            # Forward pass through the sequence
            for pos in range(len(seq_list) - k + 1):
                kmer = ''.join(seq_list[pos:pos+k]) # Make kmer
                # Check for errors, and if there are, run correction
                if check_kmer_error(kmer, kmer_frequency, threshold):
                    if args.verbose:
                        print(f'K-mer {kmer} is below the error threshold.')
                    window = seq_list[pos:pos+k] # Make window of kmer
                    new_window = make_window_correction(args, 
                        window, kmer_frequency, threshold, forward=True)
                    seq_list[pos:pos+k] = list(new_window) # Replace old section of sequence with corrected section
                    if args.verbose:
                        print(f'Corrected nucleotides: {new_window}')
                else: # Checks for first good kmer after initial error
                    if error_pos == -1:
                        error_pos = pos # remember position of known good kmer

            # Backward pass if the first k-mer has an error; similar to forward direction
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

            # Correct the sequence and then add to copy of record
            new_sequence = ''.join(seq_list)
            corrected_record = record
            corrected_record.seq = Seq(new_sequence)
            # Write corrected sequences to output file
            output_name = 'corrected_output.' + file_format
            with open(output_name, 'a') as output_file:
                SeqIO.write(corrected_record, output_file, file_format)


def main():
    "Main function"

    # Variable for performing testing on system arguments
    test_args = None
    args = parse_arguments(test_args)
    if args.verbose:
        print(args)
    
    # Open the corrected file in overwrite mode
    corrected_file_name = "corrected_output."+args.format
    with open(corrected_file_name, 'w'):
        pass
    
    # Get Kmer frequency, plot it, and then check for errors
    kmer_frequency = read_illumina_file(args, args.file)
    plot_kmer_frequencies(kmer_frequency, 'K-mer Frequencies Before Correction')
    check_sequences(args, kmer_frequency)

    # Calculate k-mer frequency from corrected file and re-plot
    corrected_kmer_frequency = read_illumina_file(args, corrected_file_name)
    plot_kmer_frequencies(corrected_kmer_frequency, 'K-mer Frequencies After Correction')


main()
