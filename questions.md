# Questions to ask Travis

## Dealing with Future Errors

- When dealing with a change, that causes more errors down the line, what do you do?
- Right now I am having issues where when I change a nucleotide within an erroneous kmer, that change will cause a future kmer to be a kmer that was not originally hashed. 
- This obviously causes an error, so how do I deal with it.

### Check every future kmer on change

- A possible solution that I've thought of is checking every possible future kmer to ensure it exists when changing a base

#### Example

```
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