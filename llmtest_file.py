from Bio import SeqIO
from collections import Counter, defaultdict
import argparse

def get_kmers(sequence, k):
    """Generate k-mers from a sequence."""
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]

def count_kmers(sequences, k):
    """Count k-mers across all sequences."""
    kmer_counts = Counter()
    for seq in sequences:
        kmers = get_kmers(str(seq.seq), k)
        kmer_counts.update(kmers)
    return kmer_counts

def find_correction(kmer, kmer_counts, k):
    """Find the best correction for an erroneous k-mer."""
    # If the k-mer exists, it's not an error
    if kmer_counts[kmer] > 0:
        return kmer
    
    # Find similar k-mers (1 substitution away)
    candidates = []
    bases = ['A', 'T', 'C', 'G']
    for i in range(len(kmer)):
        for base in bases:
            if kmer[i] != base:
                candidate = kmer[:i] + base + kmer[i+1:]
                if kmer_counts[candidate] > 0:
                    candidates.append((candidate, kmer_counts[candidate]))
    
    # Return the candidate with the highest count, if any
    if candidates:
        return max(candidates, key=lambda x: x[1])[0]
    
    # If no suitable correction found, return the original k-mer
    return kmer

def correct_sequence(sequence, kmer_counts, k):
    """Correct a sequence based on k-mer counts."""
    kmers = get_kmers(str(sequence.seq), k)
    corrected_seq = list(sequence.seq)
    
    for i, kmer in enumerate(kmers):
        corrected_kmer = find_correction(kmer, kmer_counts, k)
        if corrected_kmer != kmer:
            # Replace the erroneous base with the corrected one
            corrected_seq[i:i+k] = corrected_kmer
    
    sequence.seq = "".join(corrected_seq)
    return sequence

def error_correct(input_file, output_file, k=31, threshold=2):
    """Error correct an Illumina FASTA/FASTQ file using k-mers."""
    # Step 1: Read sequences from the input file
    sequences = list(SeqIO.parse(input_file, "fastq" if input_file.endswith(".fastq") else "fasta"))
    
    # Step 2: Count k-mers across all sequences
    kmer_counts = count_kmers(sequences, k)
    
    # Step 3: Correct sequences
    corrected_sequences = []
    for seq in sequences:
        corrected_seq = correct_sequence(seq, kmer_counts, k)
        corrected_sequences.append(corrected_seq)
    
    # Step 4: Write the corrected sequences to the output file
    with open(output_file, "w") as out_f:
        SeqIO.write(corrected_sequences, out_f, "fastq" if input_file.endswith(".fastq") else "fasta")

print(__name__)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Error correct Illumina sequences using k-mers")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA/FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA/FASTQ file")
    parser.add_argument("-k", "--kmer", type=int, default=31, help="Length of k-mers (default: 31)")
    parser.add_argument("-t", "--threshold", type=int, default=2, help="Threshold for k-mer frequency (default: 2)")
    
    args = parser.parse_args()
    
    # Error correct the input file
    error_correct(args.input, args.output, k=args.kmer, threshold=args.threshold)
