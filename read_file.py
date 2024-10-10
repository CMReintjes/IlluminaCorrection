from msilib import sequence
import Bio
from Bio.Seq import Seq
from Bio import SeqIO


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


