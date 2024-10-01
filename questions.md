# Questions to ask Travis

## Dealing with Future Errors

- When dealing with a change, that causes more errors down the line, what do you do?
- Right now I am having issues where when I change a nucleotide within an erroneous kmer, that change will cause a future kmer to be a kmer that was not originally hashed.\
- This obviously causes an error, so how do I deal with it. \
- When a kmer is below the threshold, AND one of its variants has a higher frequency, do you HAVE to change it.

### Example

```txt
                                                               TGGCAAACTGCCCGGAAAGAA
                                                               TGGCAAACTGCCCGGAAAGAC
CAAGACAGTGCCGAGCACCTAGAAGACAGGGCTGCTGGAAGTGTGGCAAGCCAGGACACATCATGGCAAACTGCCCGGAAAGACTGGCAGGTTTTTGAGG
Kmer: GGCAAACTGCCCGGAAAGACT does not exist within the sequence
Kmer: GCAAACTGCCCGGAAAGACTG does not exist within the sequence
Kmer: CAAACTGCCCGGAAAGACTGG does not exist within the sequence
Kmer: AAACTGCCCGGAAAGACTGGC does not exist within the sequence
Kmer: AACTGCCCGGAAAGACTGGCA does not exist within the sequence
Kmer: ACTGCCCGGAAAGACTGGCAG does not exist within the sequence
Kmer: CTGCCCGGAAAGACTGGCAGG does not exist within the sequence
Kmer: TGCCCGGAAAGACTGGCAGGT does not exist within the sequence
Kmer: GCCCGGAAAGACTGGCAGGTT does not exist within the sequence
Kmer: CCCGGAAAGACTGGCAGGTTT does not exist within the sequence
Kmer: CCGGAAAGACTGGCAGGTTTT does not exist within the sequence
Kmer: CGGAAAGACTGGCAGGTTTTT does not exist within the sequence
Kmer: GGAAAGACTGGCAGGTTTTTG does not exist within the sequence
Kmer: GAAAGACTGGCAGGTTTTTGA does not exist within the sequence
Kmer: AAAGACTGGCAGGTTTTTGAG does not exist within the sequence
Kmer: AAGACTGGCAGGTTTTTGAGG does not exist within the sequence
```

### Check every future kmer on change

- A possible solution that I've thought of is checking every possible future kmer to ensure it exists when changing a base.

## Using Hamming Distance for Correction

- Instead of changing only the last base, I could using hamming distance to determine the fewest number of changes anywhere in a kmer that need to be made to correct the kmer.

- This prevents from having to deal with the edge-case of errors in the initial kmer.
