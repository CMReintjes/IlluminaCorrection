# Process Flow

In this assignment, the task is to develop a basic Illumina read error correction software that leverages _k-mer_ frequencies to identify and fix errors in sequencing reads. Here’s how you can break down the error correction process into a series of steps:

## 1. **Calculate k-mer Frequencies**

- **Input:** A set of DNA reads.
- **Output:** A hash map or dictionary that stores the counts of each _k-mer_ across all the reads.
- For each read, generate all possible _k-mers_ (substrings of length `k`). Use a hash map to count the occurrences of each _k-mer_.

```python
def calculate_kmer_frequencies(reads, k):
    kmer_counts = {}
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer in kmer_counts:
                kmer_counts[kmer] += 1
            else:
                kmer_counts[kmer] = 1
    return kmer_counts
```

## 2. **Determine the Error Threshold**

- Since the genome coverage is `X = 30`, the average _k-mer_ count for true sequences should be around 30.
- The threshold can be set to a low count (e.g., 1 or 2) to identify _k-mers_ that likely contain errors. This means any _k-mer_ with a count below this threshold is considered an error.

## 3. **Identify Potential Errors in the Reads**

- For each read, slide through it and generate _k-mers_. If a _k-mer_ has a count below the error threshold, the base in that region is flagged as a potential error.

## 4. **Correct Errors**

- For each flagged error position, change the base to one of the three other possible nucleotides (A, T, G, C).
- Generate new _k-mers_ using the modified base and check their counts in the _k-mer_ frequency hash map.
- Select the base that results in _k-mers_ with the highest counts.

```python
def correct_errors(reads, k, kmer_counts, error_threshold=2):
    corrected_reads = []

    for read in reads:
        read_list = list(read)  # Convert to list for easy modification
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            if kmer_counts.get(kmer, 0) <= error_threshold:
                # Possible error detected; find the base in the k-mer causing the issue
                for j in range(k):  # Loop through each base in the k-mer
                    original_base = read[i+j]
                    max_count = 0
                    best_base = original_base

                    # Try replacing with each possible base
                    for base in 'ATGC':
                        if base != original_base:
                            # Create a new k-mer with the modified base
                            modified_kmer = list(kmer)
                            modified_kmer[j] = base
                            modified_kmer_str = ''.join(modified_kmer)

                            # Check the count of the new k-mer
                            count = kmer_counts.get(modified_kmer_str, 0)
                            if count > max_count:
                                max_count = count
                                best_base = base

                    # Replace the original base with the best base if it's different
                    if best_base != original_base:
                        read_list[i+j] = best_base

        corrected_reads.append(''.join(read_list))

    return corrected_reads
```

## 5. **Run the Correction Process**

- Load your reads and select an appropriate _k_ (e.g., `k = 21`).
- Calculate _k-mer_ frequencies for all reads.
- Iterate through each read to identify and correct errors based on _k-mer_ counts.

```python
# Sample usage
reads = ["GAGCGAAAGGGCAGACGAAGAAGGACTCCAAGGAAAACTAAGACTCCTTCGCCTTCTGCACCAGACAAGTGAGTATGGAGCCTGGTAGGAATCAGCTGTT"]
k = 21

kmer_counts = calculate_kmer_frequencies(reads, k)
corrected_reads = correct_errors(reads, k, kmer_counts)

print(corrected_reads)
```

## Summary of the Steps

1. **Calculate _k-mer_ frequencies** for the given set of reads.
2. **Set an error threshold** (e.g., 1 or 2) to identify error-prone _k-mers_.
3. **Flag errors** in reads by checking _k-mer_ counts against the threshold.
4. **Correct errors** by replacing the suspect base with the one that results in the highest _k-mer_ counts.
5. **Return corrected reads**.

## Notes

- Adjusting parameters like _k_ and the error threshold can significantly affect the performance and accuracy of the correction.
- The algorithm works well for small reads and error rates typical of Illumina data, especially when sequencing depth is high (coverage `X ≥ 30`).
