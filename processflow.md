# Process Flow

In this assignment, the task is to develop a basic Illumina read error correction software that leverages _k-mer_ frequencies to identify and fix errors in sequencing reads. Here’s how you can break down the error correction process into a series of steps:

## 1. **Calculate k-mer Frequencies**

- **Input:** A set of DNA reads.
- **Output:** A hash map or dictionary that stores the counts of each _k-mer_ across all the reads.
- For each read, generate all possible _k-mers_ (substrings of length `k`). Use a hash map to count the occurrences of each _k-mer_.

## 2. **Determine the Error Threshold**

- Since the genome coverage is `X = 30`, the average _k-mer_ count for true sequences should be around 30.
- The threshold can be set to a low count (e.g., 1 or 2) to identify _k-mers_ that likely contain errors. This means any _k-mer_ with a count below this threshold is considered an error.

## 3. **Identify Potential Errors in the Reads**

- For each read, slide through it and generate _k-mers_. If a _k-mer_ has a count below the error threshold, the base in that region is flagged as a potential error.

## 4. **Correct Errors**

- For each flagged error position, change the base to one of the three other possible nucleotides (A, T, G, C).
- Generate new _k-mers_ using the modified base and check their counts in the _k-mer_ frequency hash map.
- Select the base that results in _k-mers_ with the highest counts.

## 5. **Run the Correction Process**

- Load your reads and select an appropriate _k_ (e.g., `k = 21`).
- Calculate _k-mer_ frequencies for all reads.
- Iterate through each read to identify and correct errors based on _k-mer_ counts.

## Summary of the Steps

1. **Calculate _k-mer_ frequencies** for the given set of reads.
2. **Set an error threshold** (e.g., 1 or 2) to identify error-prone _k-mers_.
3. **Flag errors** in reads by checking _k-mer_ counts against the threshold.
4. **Correct errors** by replacing the suspect base with the one that results in the highest _k-mer_ counts.
5. **Return corrected reads**.

## Notes

- Adjusting parameters like _k_ and the error threshold can significantly affect the performance and accuracy of the correction.
- The algorithm works well for small reads and error rates typical of Illumina data, especially when sequencing depth is high (coverage `X ≥ 30`).
