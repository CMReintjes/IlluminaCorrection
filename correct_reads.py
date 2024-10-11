from Bio.Seq import Seq
from Bio import SeqIO
import matplotlib.pyplot as plt
import argparse


def parse_arguments(test_args):
    "Parse command line arguments as input"
    parser = argparse.ArgumentParser(description='Error correction for Illumina sequencing data using k-mer counts')

    # Positional arguments
    parser.add_argument('file', type=str, help='Illumina reads file for correction')
    parser.add_argument('format', type=str.lower, choices=['fasta', 'fastq'], help='Format of the Illumina file (.fasta or .fastq)')
    parser.add_argument('kmer_length', type=int, help='Length of the k-mer')
    parser.add_argument('threshold', type=int, help='Error threshold for correction')

    # Optional flag arguments
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose output')
    parser.add_argument('-s', '--save-plot', action='store_true')

    if test_args:
        return parser.parse_args(test_args)
    else:
        return parser.parse_args()


def plot_kmer_frequencies(kmer_frequency, title, close=False):
    "Plots a histogram of k-mer frequencies."
    plt.figure(figsize=(10, 6))
    plt.hist(kmer_frequency.values(), bins=50)
    plt.xlabel('K-mer Frequency')
    plt.ylabel('Count')
    plt.title(title)
    if close:
        plt.savefig(title+".png", dpi=300)
        plt.close()
    else:
        plt.show()


def get_kmer_frequency(kmer, frequency):
    "Checks if kmer exists within the hash, and updates count"
    if kmer in frequency: # If key exists add one to frequency
        frequency[kmer] += 1
    else: # If new kmer then add with frequency 1
        frequency[kmer] = 1
    return


def read_illumina(filename, format, k, verbose=False):
    """
    Read and parse the illumina file and return kmer hash

    Args:
        filename: name of the input file
        format: file format for biopython parser
        k: length of the kmers

    Returns:
        kmer_frequency: hash of kmers and count frequency
    """
    kmer_frequency = {}

    if verbose:
        print(f'Opening file: {filename}')
    # Open Illumina File
    with open(filename, 'r') as input_file:
        if verbose: # Verbose output
            print(f'Parsing file as {format} format...')
        # Iterate through records
        for record in SeqIO.parse(input_file, format):
            sequence = str(record.seq)
            if verbose: # Verbose output
                print(f'Processing record: {record.id}')
            # Get kmer frequency using dictionary
            for pos in range(len(sequence)-k+1):
                get_kmer_frequency(sequence[pos:pos+k], kmer_frequency)
    return kmer_frequency


def check_kmer_error(kmer, kmer_frequency, threshold):
    "Checks if a k-mer is below the error threshold. Returns True if it has an error."
    try: # Check if kmer is below the error threshold
        if kmer_frequency[kmer] < threshold:
            return True
        else:
            return False
    except KeyError: # Error handling just in case
        print(f"Error encountered with kmer: {kmer}. Key missing from hash.")
        return True
    

def fix_window(window, k, threshold, frequency, reverse=False):
    """
    Use the sliding window to check the kmer for possible base changes, and check validity of future affected kmers.

    Args:
        window: list of nucleotides to correct (len=2k-1)
        k: kmer length
        threshold: error threshold
        frequency: dictionary of kmer frequency
        reverse: boolean for changing to reverse check direction

    Returns:
        corrected: corrected kmer (len=k)
    """
    corrected = list(window)
    all_kmers = True
    if reverse: # Set values for reverse checking sequences
        pos = len(window)-k
        start_pos = len(window)-k
        end_pos = -1
        step = -1
        original = ''.join(window[start_pos:start_pos+k+1])
    else: # Default values for forward direction
        pos = k-1
        start_pos = 0
        end_pos = len(window)-k+1
        step = 1
        original = ''.join(window[start_pos:start_pos+k])
    
    max_freq = frequency[original]
    bases = list("ATGC")
    for b in bases: # Iterate through base options
        all_kmers = True
        window[pos] = b # Change the base at error location
        # Iterate through future kmers to check validity
        for n in range(start_pos, end_pos, step):
            kmer = "".join(window[n:n+k])
            # If kmer doesn't exist break out
            if kmer not in frequency:
                all_kmers = False
                break

        # Make kmer with changed original
        kmer = list(original)
        kmer[pos] = b
        kmer = ''.join(kmer)
        # Check that kmer is above the frequency of the last and over threshold
        if all_kmers and kmer in frequency:
            if frequency[kmer] > max_freq and frequency[kmer] > threshold:
                corrected[pos] = b

    return corrected[start_pos:start_pos+k]


def check_sequences(args, kmer_frequency):
    """
    Parses an Illumina sequencing file, and then corrects erroneous bases.
    
    First parses the file and iterates through sequences. Sequences are made into a list
    and then split into kmers. If kmer is erroneous, sliding window is then made of 
    length 2k-1 to contain all kmers that would be affected by the change. If there 
    is an error within the first kmer, they are skipped until the first valid kmer, 
    and are then checked at the end using a reverse pass through the initial area.

    Args:
        args: argparse arguments object
        kmer_frequency: dictionary of kmer frequencies
    
    Outputs:
        Corrected sequences to corrected file
    """
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
                        print(f'K-mer {kmer} is below the error threshold: {kmer_frequency[kmer]}.')
                    window = list(seq_list[pos:pos+(2*k)-1]) # Make window of kmer
                    new_window = fix_window(list(window), k, threshold, kmer_frequency)
                    seq_list[pos:pos+k] = list(new_window) # Replace old section of sequence with corrected section
                    if args.verbose:
                        print(f'Corrected kmer: {''.join(new_window)}')
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
                        adj_pos = max(0, pos-k)
                        end_pos = pos+(k)
                        window = list(seq_list[adj_pos:end_pos])
                        new_window = fix_window(list(window), k, threshold, kmer_frequency, reverse=True)
                        seq_list[pos:pos+k] = list(new_window)
                        if args.verbose:
                            print(f'Backward pass corrected kmer: {''.join(new_window)}')

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
    kmer_frequency = read_illumina(args.file, args.format, args.kmer_length)
    plot_kmer_frequencies(kmer_frequency, 'K-mer Frequencies Before Correction', args.save_plot)
    check_sequences(args, kmer_frequency)

    # Calculate k-mer frequency from corrected file and re-plot
    corrected_kmer_frequency = read_illumina(corrected_file_name, args.format, args.kmer_length)
    plot_kmer_frequencies(corrected_kmer_frequency, 'K-mer Frequencies After Correction', args.save_plot)


main()
