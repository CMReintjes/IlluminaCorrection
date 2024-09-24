# Ideas

## Kmer Problem

### Kmer List

- When you have a kmer with an error, grab all future kmers including the error nucleotide.
- Do this in a list (length is kmer length).
- Then change the nucleotide in all of them simultaneously and check which have the best overall score.

- A possibility (to check with TB) is to take an average of the frequency across this list and use the highest average

### Iterative Kmer Correction

- Once kmers have been changed, and subsequent kmers don't exist within the dictionary, rehash the kmers.
- This can be combined with the previous idea.

#### Per Sequence Kmer Correction

- The previously mentioned idea can be done on a local basis on just the current sequence.
- Iterate through the sequence again, subtracting the counts of kmers that exist within the original (un-corrected sequence) and then rehash just that one sequence.
- This can be done after every change.
