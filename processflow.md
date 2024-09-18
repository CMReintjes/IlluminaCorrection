# Process Flow

## 1. **Calculate k-mers and their Frequency (k-counts)**

- The DNA sequence is split into small overlapping pieces called _k-mers_ of length "k."
- Each _k-mer_ has a count (or frequency), called the k-count, representing how often that _k-mer_ appears in the dataset.
- If a _k-mer_ has a **high k-count**, it’s considered a "solid" (likely correct) _k-mer_. If it has a **low k-count**, it’s considered "weak" (likely to have errors).
- Because sequencing data can be huge, a hashed map (a fast way of storing data) is used to store each _k-mer_ and its positions in the read.

To store _k-mer_ counts efficiently, a **hash map (hash table)** can be used. A hash map allows quick lookup, insertion, and counting of _k-mers_ by using the _k-mer_ sequence as the key. Here's how you can implement a hash map for _k-mer_ counts step-by-step:

### 1. **Initialize the Hash Map**

- Create an empty hash map (often implemented as a dictionary in Python, a `HashMap` in Java, or an `unordered_map` in C++).
- In this hash map, the keys will be the _k-mers_ (strings of DNA bases of length `k`), and the values will be their counts (integers).

### 2. **Iterate Over Each DNA Read**

- For each DNA read in your dataset, extract all possible _k-mers_ of length `k`.
- In a read of length `L`, you can generate `L - k + 1` _k-mers_. For example, for a sequence `ACGTG` and `k = 3`, the _k-mers_ are `ACG`, `CGT`, and `GTG`.

### 3. **Hash Each k-mer and Update Its Count**

- For each _k-mer_, check if it already exists in the hash map:
  - If the _k-mer_ is **not in the hash map**, add it to the hash map with an initial count of 1.
  - If the _k-mer_ **exists in the hash map**, increment its count by 1.
- Here's a Python example for clarity:

```python
def count_kmers(sequence, k):
    kmer_counts = {}  # Initialize an empty hash map (dictionary)

    # Iterate through the sequence to extract all k-mers
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]  # Get the k-mer of length k starting at index i

        # Update the k-mer count in the hash map
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
        else:
            kmer_counts[kmer] = 1

    return kmer_counts
```

- This function takes a DNA sequence and a _k_-value, iterates through the sequence to generate _k-mers_, and updates their counts in the hash map.

### 4. **Handle Large Data (Memory Management)**

- For large datasets, you may need to use more memory-efficient methods:
  - **Sliding Window**: Process the reads in chunks rather than loading the entire dataset into memory.
  - **Disk-Based Storage**: Use disk-based key-value stores (like LevelDB) if the in-memory hash map becomes too large.
  - **Hash Collisions**: Make sure to use a good hash function (most built-in hash functions are suitable) to minimize collisions and optimize access times.

### 5. **Using Hashing for Fast Lookup**

- The hash function takes each _k-mer_ (string) and computes an index to store it in the hash map, allowing quick insertion and retrieval.
- A well-designed hash function ensures that each _k-mer_ is distributed evenly across the hash map, making counting fast and efficient.

### Kmer Frequency Summary

- **Keys:** The _k-mers_ (substrings of length `k`).
- **Values:** The counts of each _k-mer_.
- Iterate over the DNA sequence to extract _k-mers_ and update the counts in the hash map.
- This approach allows efficient counting and storage of _k-mers_, making it a practical solution for error correction in large sequencing datasets.

## 2. **Establish the Error Threshold**

- This step separates solid _k-mers_ (correct) from weak _k-mers_ (incorrect).
- The error threshold is a k-count value that distinguishes between solid and weak _k-mers_.
- To find this threshold, the frequency of each k-count is examined. Typically, the first minimum in the frequency distribution of k-counts is used to set the threshold.
- In simpler terms, this means finding the point where the count of _k-mers_ shifts from being common (solid) to rare (weak).
- If there is a long section in the frequency distribution with consecutive zeros (no _k-mers_), this is used to set the error threshold.

To establish the error threshold, the method involves analyzing the frequency distribution of _k-mer_ counts to distinguish between solid (correct) and weak (error-prone) _k-mers_. Here’s how it is done step-by-step:

### 1. **Calculate the Frequency Distribution of k-mer Counts**

- After counting how many times each _k-mer_ appears in the sequencing data, you create a distribution (a table or graph) that shows how many _k-mers_ exist for each possible count.
- For example, you might find many _k-mers_ that appear 20 times, fewer that appear 5 times, and very few that appear just once or twice.

### 2. **Identify the Minimum Point in the Frequency Distribution**

- The frequency distribution usually has two distinct parts:
  1.  **High-frequency region:** Representing _k-mers_ that appear often (solid, correct _k-mers_).
  2.  **Low-frequency region:** Representing _k-mers_ that appear infrequently (likely to be errors).
- These two regions are typically separated by a "valley" or a dip in the frequency distribution. This valley indicates a point where the frequency of _k-mers_ transitions from being common (solid) to rare (weak).
- The first minimum in this valley is identified as a potential error threshold.

### 3. **Handle Irregularities in Distribution (e.g., Consecutive Zeros)**

- In some cases, there might be sections in the distribution where the count drops to zero for consecutive values. This happens when there are _k-count_ values with no corresponding _k-mers_.
- The threshold is often set at the end of the first suitably long section of consecutive zeros. This section marks a clear boundary between common and rare _k-mers_.
- The length of this zero section helps define the error threshold, as it signals that _k-mers_ beyond this point are too rare to be correct.

### 4. **Define the Error Threshold**

- The error threshold (`ter`) is defined as the _k-mer_ count value at the end of this zero section or the first significant minimum point in the distribution.
- This threshold is then used to classify _k-mers_: those with counts below the threshold are considered weak (potential errors), while those with counts equal to or above the threshold are treated as solid (likely correct).

### Error Threshold Summary

- The error threshold is found by analyzing the distribution of _k-mer_ counts.
- It is set at the point where the distribution has its first minimum or a prolonged drop (e.g., consecutive zeros).
- This threshold helps separate solid _k-mers_ (correct) from weak _k-mers_ (likely errors), which can then be used to identify and correct sequencing errors.

## 3. **Identify Error Regions in the Sequence**

- The sequence is scanned for regions where _k-mers_ have k-counts below the error threshold. These regions are flagged as potential error zones.
- The flagged regions are grouped into clusters using a technique called variable bandwidth mean-shift clustering (to identify error areas more accurately).
- After identifying these error sections, they can be extended in both directions to include adjacent _k-mers_ that might belong to the same error cluster.

## 4. **Correct Errors in the Error Sections**

- Once error regions are identified, the system attempts to correct them.
- It analyzes the surrounding solid _k-mers_ to suggest a correction, replacing the weak _k-mers_ with more likely, high-frequency ones.
- This process improves the accuracy of the DNA sequence data by fixing these error-prone regions.

## Full Summary

- **Phase 1:** Calculate how often each _k-mer_ appears in the sequence.
- **Phase 2:** Find a threshold that separates common (solid) _k-mers_ from rare (weak) ones.
- **Phase 3:** Identify regions in the DNA sequence where there are likely errors based on the threshold.
- **Phase 4:** Correct the errors using information from the solid (correct) _k-mers_ around them.

This method makes use of the natural repetition of DNA sequences to distinguish between real sequence data and sequencing errors.
