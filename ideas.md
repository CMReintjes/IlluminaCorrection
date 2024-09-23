# Ideas

## Kmer Problem

When you have a kmer with an error, grab all future kmers including the error nucleotide. Do this in a list (length is kmer length). Then change the nucleotide in all of them simultaneously and check which have the best overall score.