from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import argparse


def parse_arguments(test_args):
    '''
    Parse command line arguments as input.

    This function sets up and parses command-line arguments for the script, allowing users to specify the input Illumina file, file format, k-mer length, error threshold, and verbosity.

    Args:
        test_args (list): List of arguments for testing purposes.

    Returns:
        Namespace: Parsed arguments including file, format, kmer_length, threshold, and verbose.
    '''
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
    '''
    Parse the Illumina file and get k-mer frequency from sequences.

    This function reads an Illumina sequencing file, processes each sequence, and computes the frequency of all k-mers within the file. The k-mer frequencies are collected into a dictionary.

    Args:
        args (Namespace): Parsed command line arguments.
        file_name (str): Path to the Illumina file.

    Returns:
        dict: A dictionary of k-mer frequencies.
    '''
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
            # Get k-mer frequency using dictionary
            kmer_frequency = get_kmer_frequency(str(record.seq), k, kmer_frequency)
    return kmer_frequency


def get_kmer_frequency(sequence, k, kmer_frequency):
    '''
    Get k-mers from the sequence and update the frequency dictionary.

    This function extracts all possible k-mers from a given sequence and updates the provided dictionary with their frequencies. If the k-mer is new, it is added to the dictionary; otherwise, its frequency is incremented.

    Args:
        sequence (str): DNA sequence from which to extract k-mers.
        k (int): Length of the k-mer.
        kmer_frequency (dict): Dictionary containing k-mer frequencies.

    Returns:
        dict: Updated k-mer frequency dictionary.
    '''
    for pos in range(len(sequence)-k+1):
        # Create kmers of length k using string positions
        kmer = sequence[pos:pos+k]
        if kmer in kmer_frequency: # If key exists add one to frequency
            kmer_frequency[kmer] += 1
        else: # If new kmer then add with frequency 1
            kmer_frequency[kmer] = 1

    return kmer_frequency


def plot_kmer_frequencies(kmer_frequency, title):
    '''
    Plot a histogram of k-mer frequencies.

    This function generates a histogram that visualizes the distribution of k-mer frequencies, helping to identify common and rare k-mers in the dataset.

    Args:
        kmer_frequency (dict): Dictionary containing k-mer frequencies.
        title (str): Title for the plot.
    '''
    plt.figure(figsize=(10, 6))
    plt.hist(kmer_frequency.values(), bins=50)
    plt.xlabel('K-mer Frequency')
    plt.ylabel('K-mer Count')
    plt.title(title)
    plt.show()


def check_kmer_error(kmer, kmer_frequency, threshold):
    '''
    Check if a k-mer is below the error threshold.

    This function checks whether the frequency of a given k-mer is below the specified error threshold, indicating that the k-mer might be erroneous.

    Args:
        kmer (str): The k-mer to be checked.
        kmer_frequency (dict): Dictionary containing k-mer frequencies.
        threshold (int): Error threshold for the k-mer frequency.

    Returns:
        bool: True if the k-mer has an error (i.e., frequency is below threshold), False otherwise.
    '''
    try: # Check if kmer is below the error threshold
        if kmer_frequency[kmer] < threshold:
            return True
    except KeyError: # Error handling just in case
        print(f'Kmer: {kmer} not in hash. Expected poor previous base change')
        return True  


def make_window_correction(args, seq_list, start_pos, k, kmer_frequency, threshold, forward):
    '''
    Make corrections to an erroneous base by trying all possible base changes.

    This function attempts to correct an erroneous base by substituting it with each possible nucleotide (A, T, C, G). It checks if the corrected k-mer is present in the frequency dictionary and returns the first successful correction.

    Args:
        args (Namespace): Parsed command line arguments.
        seq_list (list): List of characters representing the sequence.
        start_pos (int): Starting position of the k-mer in the sequence.
        k (int): Length of the k-mer.
        kmer_frequency (dict): Dictionary containing k-mer frequencies.
        threshold (int): Error threshold for the k-mer frequency.
        forward (bool): Direction of the correction (True for forward, False for backward).

    Returns:
        list: Corrected sequence list.
    '''
    nucleotides = list('ATCG') # List of nucleotides for base changes
    position = start_pos + (k - 1) if forward else start_pos
    original = seq_list[position]

    corrected_seq_list = list(seq_list) # Make copy of original
    for nuc in nucleotides:
        if nuc != original:
            corrected_seq_list[position] = nuc
            all_kmers_exist = True
            
            # Check all kmers affected by this change, ensuring they exist in the dictionary
            for pos in range(start_pos, min(start_pos + k, len(seq_list) - k + 1)):
                new_kmer = ''.join(corrected_seq_list[pos:pos + k])
                if len(new_kmer) == k and new_kmer not in kmer_frequency:
                    all_kmers_exist = False
                    break

            if all_kmers_exist:
                if args.verbose:
                    print(f'Correction successful for base at position {position}: {nuc}')
                return corrected_seq_list

    # If no correction was successful, revert the change
    if args.verbose:
        print(f'Correction failed for base at position {position}')
    return seq_list


def check_sequences(args, kmer_frequency):
    '''
    Parse an Illumina sequencing file and correct erroneous bases.

    This function reads sequences from an Illumina sequencing file and corrects erroneous k-mers by performing forward and backward passes. It generates corrected sequences and writes them to an output file.

    Args:
        args (Namespace): Parsed command line arguments.
        kmer_frequency (dict): Dictionary containing k-mer frequencies.
    '''
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
                    new_seq_list = make_window_correction(args, seq_list, pos, k, kmer_frequency, threshold, forward=True)
                    seq_list[pos:pos+k] = new_seq_list[pos:pos+k] # Replace old section of sequence with corrected section
                    if args.verbose:
                        print(f'Corrected nucleotides: {"".join(new_seq_list[pos:pos+k])}')
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
                        new_seq_list = make_window_correction(args, seq_list, pos, k, kmer_frequency, threshold, forward=False)
                        seq_list[pos:pos+k] = new_seq_list[pos:pos+k]
                        if args.verbose:
                            print(f'Backward pass corrected nucleotides: {"".join(new_seq_list[pos:pos+k])}')

            # Correct the sequence and then add to copy of record
            new_sequence = ''.join(seq_list)
            corrected_record = record
            corrected_record.seq = Seq(new_sequence)
            # Write corrected sequences to output file
            output_name = 'corrected_output.' + file_format
            with open(output_name, 'a') as output_file:
                SeqIO.write(corrected_record, output_file, file_format)


def main():
    '''
    Main function that orchestrates the error correction for Illumina sequencing data.

    This function handles the entire process of reading input arguments, computing initial k-mer frequencies, correcting erroneous k-mers, and saving the corrected sequences. It also generates visualizations of k-mer frequencies before and after correction.
    '''

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
