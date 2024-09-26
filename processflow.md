# Process Flow

The objective of this assignment is to develop a basic Illumina read error correction software that uses _k-mer_ frequencies and an error threshold to identify and fix errors in sequenced reads.

## 1. **Calculate k-mer Frequencies**

- **Input:** A set of DNA reads.
- **Output:** A hash map or dictionary that stores the counts of each _k-mer_ across all the reads.
- For each read, generate all possible _k-mers_ (substrings of length `k`). Use a hash map to count the occurrences of each _k-mer_.

### 1.1. **Create Histogram of Frequencies**

- Create a histogram of the kmer frequency distribution
- Use `seaborn` or `pyplot`

## 2. **Identify Potential Errors in the Reads**

- For each read, slide through it and generate _k-mers_. If a _k-mer_ has a count below the error threshold, the base in that region is flagged as a potential error.

## 3. **Correct Errors**

- For each flagged error position, change the base to one of the three other possible nucleotides (A, T, G, C).
- Generate new _k-mers_ using the modified base and check their counts in the _k-mer_ frequency hash map.
- Select the base that results in _k-mers_ with the highest counts.

### 3.1. **All-Kmer Nucleotide Check**

- When checking the kmer, check all of the following to ensure that the change all kmers that contain that error nucleotide exist within the hash.

## 4. **Check the Correction Process**

- Load the corrected reads and re-hash the kmers.
- Check for any remaining errors within the reads.

### 4.1. File Output and Histogram

- Output the corrected reads to a file of same format as input
- Output a histogram of the kmer distribution

## Summary of the Steps

1. **Calculate _k-mer_ frequencies** for the given set of reads.
2. **Flag errors** in reads by checking _k-mer_ counts against the threshold.
3. **Correct errors** by replacing the suspect base with the one that results in the highest _k-mer_ counts.
4. **Return corrected reads**.

## Notes

- Adjusting parameters like _k_ and the error threshold can significantly affect the performance and accuracy of the correction.
- The algorithm works well for small reads and error rates typical of Illumina data, especially when sequencing depth is high (coverage `X ≥ 30`).
