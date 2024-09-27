# Process Flow

K-mer based error correction is a popular method to reduce sequencing errors in Illumina reads, which tend to have relatively low error rates compared to other technologies. Illumina reads are typically provided in FASTA or FASTQ formats. The process of K-mer based error correction can be implemented in Python, and the use of Biopython for reading and writing the data can help streamline the workflow.

Here’s a typical process flow:

## 1. **Read Input Data**

- **Input formats**: The Illumina reads are provided in **FASTA** or **FASTQ** formats.
- **Biopython** can be used to parse the files:
  - For **FASTQ**: `Bio.SeqIO.parse("input.fastq", "fastq")`
  - For **FASTA**: `Bio.SeqIO.parse("input.fasta", "fasta")`

   This step allows you to read sequences into memory to perform K-mer analysis.

## 2. **Build K-mer Frequency Table**

- **K-mers**: These are all substrings of length *k* found within each sequence. You need to determine an appropriate value for *k* (e.g., 31, 51) based on the read length and the complexity of the genome.
- Iterate over the sequences to extract all possible k-mers and maintain a count of how frequently each k-mer occurs in the dataset.

   Example:

   ```python
   from collections import defaultdict
   
   k = 31  # Example value for k
   kmer_counts = defaultdict(int)
   
   for record in Bio.SeqIO.parse("input.fastq", "fastq"):
       seq = str(record.seq)
       for i in range(len(seq) - k + 1):
           kmer = seq[i:i+k]
           kmer_counts[kmer] += 1
   ```

## 3. **Identify Likely Error K-mers**

- **Thresholding**: Typically, high-frequency k-mers represent true biological sequences, while low-frequency k-mers are likely to contain sequencing errors.
- **Cutoff**: Establish a frequency cutoff to classify k-mers as erroneous or correct. For example, k-mers that occur fewer than 2–5 times may be considered errors.
- **Set of trusted k-mers**: Create a set of "trusted" k-mers based on the frequency threshold.

Determining **trusted k-mers** versus **erroneous k-mers** (based on a cutoff frequency) is a critical step in k-mer-based error correction. This step involves setting a **cutoff threshold** for k-mer frequency, which distinguishes between k-mers likely originating from true biological sequences and those resulting from sequencing errors.

### **Process for Determining Trusted K-mers vs. Cutoff:**

1. **K-mer Frequency Distribution**:
   - After extracting k-mers from your sequencing reads, count how often each k-mer appears.
   - K-mers from the correct sequences (i.e., non-erroneous reads) tend to occur more frequently, while erroneous k-mers (which arise from random sequencing errors) typically appear at low frequencies.

   You can represent the k-mer frequencies as a histogram, where:
   - The **x-axis** shows the k-mer frequency (how many times each k-mer appears).
   - The **y-axis** shows the number of distinct k-mers with that frequency.

2. **Identifying Error-Generated K-mers**:
   - **Low-frequency k-mers** are likely to be errors, since sequencing errors are random and do not repeat often in the reads.
   - Typically, k-mers that occur only once or a few times (depending on the dataset size and sequencing depth) are suspected to be errors.

3. **Determine the Cutoff**:
   - You need to establish a **cutoff frequency** to separate trusted k-mers from erroneous k-mers. This is often done by examining the frequency distribution:
     - **Error-prone k-mers**: K-mers with very low frequencies (below the cutoff) are considered erroneous.
     - **Trusted k-mers**: K-mers with frequencies above the cutoff are treated as trusted, meaning they are likely correct.

   A common strategy is to set a **cutoff frequency** based on the characteristics of the dataset, e.g., the depth of sequencing and the organism’s genome size.

4. **Setting the Cutoff Frequency**:
   - **Visual Approach (Using Histogram)**: Plot a histogram of k-mer frequencies. The frequency distribution often shows a **bimodal pattern**:
     - The **left peak** represents low-frequency k-mers, which are mostly sequencing errors.
     - The **right peak** represents high-frequency k-mers, which are likely true biological sequences.
   - Choose a **cutoff** at the dip between the two peaks to divide the low-frequency error k-mers from the trusted ones. The exact position of this cutoff depends on the specific data and noise in the sequencing.

   Example:
   - K-mers with frequency less than **2** or **5** might be considered errors in a typical Illumina dataset, but this can vary.

5. **Automated Methods for Setting the Cutoff**:
   - **Fixed Threshold**: Use a predefined cutoff (e.g., 3 or 5) based on prior knowledge of the sequencing technology.
   - **Dynamic Threshold (Automatic Estimation)**: Tools like **KMC** and **Jellyfish** that perform k-mer counting may automatically estimate a cutoff based on the k-mer frequency distribution.
     - These tools typically assume that true k-mers occur at higher frequencies and will suggest a cutoff that separates error-generated k-mers from correct ones.

6. **Additional Considerations**:
   - **Genome Coverage**: High sequencing coverage means even low-frequency k-mers could be correct. In such cases, a higher cutoff might be necessary.
   - **Sequencing Depth**: For shallow sequencing, a lower cutoff may be more appropriate to retain more correct k-mers.
   - **Genome Complexity**: For organisms with high repeat content or large genomes, more k-mers may overlap, requiring a more careful selection of the cutoff.

### **Example Workflow for Determining Trusted K-mers**

1. **Count K-mer Frequencies**:

   ```python
   from collections import defaultdict

   kmer_counts = defaultdict(int)
   k = 31  # Example k-mer size

   for record in Bio.SeqIO.parse("input.fastq", "fastq"):
       seq = str(record.seq)
       for i in range(len(seq) - k + 1):
           kmer = seq[i:i+k]
           kmer_counts[kmer] += 1
   ```

2. **Plot K-mer Frequency Distribution**:
   Use `matplotlib` to visualize the frequency distribution of k-mers to guide the cutoff selection.

   ```python
   import matplotlib.pyplot as plt

   # Extract k-mer frequencies
   frequencies = list(kmer_counts.values())

   # Plot the histogram of k-mer frequencies
   plt.hist(frequencies, bins=100, log=True)  # Log scale to emphasize low frequencies
   plt.xlabel("K-mer Frequency")
   plt.ylabel("Number of K-mers")
   plt.show()
   ```

3. **Set a Cutoff**:
   - After examining the histogram, choose a cutoff point that separates erroneous k-mers (low frequency) from trusted k-mers (high frequency).
   - For example, set `cutoff = 5` to discard all k-mers with a frequency of less than 5:

     ```python
     trusted_kmers = {kmer for kmer, count in kmer_counts.items() if count >= 5}
     ```

4. **Use Trusted K-mers for Correction**:
   - Now, only the k-mers with a frequency greater than or equal to the cutoff are considered trusted.
   - Use these trusted k-mers for correcting erroneous reads as previously described.

### **Summary**

- **Trusted k-mers** are k-mers that appear frequently in the sequencing data and likely represent true biological sequences.
- **Erroneous k-mers** are rare and result from sequencing errors.
- The **cutoff frequency** is chosen to separate trusted k-mers from erroneous k-mers. This can be based on visual inspection of the k-mer frequency distribution, or by using automated tools.

## 4. **Error Detection**

- For each read, check if the k-mers belong to the set of trusted k-mers. If one or more k-mers in the read are not in the trusted set, the read may contain an error.

## 5. **Error Correction**

- **Correcting reads**: When an erroneous k-mer is detected, you can attempt to correct it by modifying one or more nucleotides within the k-mer to match a trusted k-mer. Common approaches include:
  - **Neighboring k-mers**: Identify correct k-mers that overlap with the erroneous k-mers and adjust the sequence.
  - **Hamming Distance**: Replace erroneous k-mers with those from the trusted set that are at a small Hamming distance (e.g., 1-2 mismatches).
- Example correction algorithm:

     ```python
     def correct_kmer(kmer, trusted_kmers):
         # Find the closest k-mer in the trusted set using Hamming distance
         closest_kmer = min(trusted_kmers, key=lambda tk: hamming_distance(kmer, tk))
         return closest_kmer
     
     def hamming_distance(s1, s2):
         return sum(c1 != c2 for c1, c2 in zip(s1, s2))
     ```

### 5.1 **Hamming Distance**

The **Hamming distance** is a metric used to measure the difference between two strings of the same length. Specifically, it counts the number of positions at which the corresponding characters (or nucleotides, in the case of DNA sequences) are different. This is useful in bioinformatics for error correction because it quantifies how dissimilar two sequences are.

#### **Example**

Consider two k-mers (short sequences) of length 6:

- `kmer1 = "ACGCTA"`
- `kmer2 = "ACGCAA"`

The **Hamming distance** between `kmer1` and `kmer2` is 1 because only the last nucleotide is different (`T` vs. `A`).

#### **Formula**

For two strings `s1` and `s2` of the same length:
\[ \text{Hamming Distance} = \sum_{i=1}^{n} \text{1 if } s1[i] \neq s2[i], \text{ else 0} \]
Where \( n \) is the length of the string.

#### **Hamming Distance in Error Correction**

- **Trusted vs. Erroneous K-mers**: In error correction, you can use Hamming distance to identify small errors in k-mers. If an erroneous k-mer differs from a trusted k-mer by a small Hamming distance (say 1 or 2 positions), it’s likely that the k-mer was affected by a sequencing error, and correcting those positions will result in a valid k-mer.
  
  For example, if a read contains the erroneous k-mer `"ACGCTA"` but the trusted set contains the k-mer `"ACGCAA"`, you can correct the error by changing the final nucleotide from `T` to `A`.

### 5.2 **Neighboring K-mers for Correction**

In a DNA read, k-mers overlap with one another, meaning consecutive k-mers share many of the same nucleotides. This property can help in error correction because the correct k-mer can often be inferred by looking at neighboring k-mers that are trusted.

#### **Concept**

- **Overlapping K-mers**: If you have a read with k-mers of length `k`, then the k-mers in the read will overlap by `k-1` bases. For example, in a sequence `"ACGCTA"`, two consecutive 4-mers are:
  - `"ACGC"` (positions 1–4)
  - `"CGCT"` (positions 2–5)

  The two k-mers share 3 bases (`"CGC"`).

#### **Neighboring in Error Correction**

1. **Identify Erroneous K-mers**: If a k-mer is not in the trusted k-mer set (low frequency), but neighboring k-mers in the read are trusted, you can attempt to correct the erroneous k-mer by making it consistent with its neighboring k-mers.

2. **Correction Example**:
   - Suppose a read contains the sequence `"ACGCTCGA"`, and the k-mers are `"ACGCT"` and `"CGCTC"`. If `"CGCTC"` is an erroneous k-mer (not found in the trusted set) but `"ACGCT"` is trusted, we could try modifying the middle nucleotide of `"CGCTC"` to match the shared part of the trusted neighboring k-mers.
   - If we find that a trusted k-mer close to `"CGCTC"` is `"ACGTC"`, and it overlaps with the previous k-mer (`"ACGCT"`), we might change the error at the last position of `"CGCTC"` to match the trusted sequence.

#### **Correction Strategy Using Neighboring K-mers**

1. **Find trusted k-mers**: Identify regions where most of the overlapping k-mers are found in the trusted set.
2. **Align erroneous k-mers**: For the erroneous k-mers, try modifying them to match neighboring trusted k-mers. This ensures consistency and accuracy.
3. **Reconstruct the read**: After making the necessary corrections to the erroneous k-mers, reconstruct the entire sequence.

#### **Example Python Code**

This simplified example shows how you might use Hamming distance and neighboring k-mers for correction:

```python
def hamming_distance(s1, s2):
    """Calculate the Hamming distance between two strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def correct_kmer(erroneous_kmer, trusted_kmers):
    """Correct an erroneous k-mer by finding the closest trusted k-mer."""
    min_distance = float('inf')
    closest_kmer = erroneous_kmer
    for trusted_kmer in trusted_kmers:
        distance = hamming_distance(erroneous_kmer, trusted_kmer)
        if distance < min_distance:
            min_distance = distance
            closest_kmer = trusted_kmer
    return closest_kmer

def correct_read_using_neighbors(read, k, trusted_kmers):
    """Attempt to correct a read by using overlapping trusted k-mers."""
    corrected_read = list(read)
    for i in range(len(read) - k + 1):
        kmer = read[i:i+k]
        if kmer not in trusted_kmers:
            corrected_kmer = correct_kmer(kmer, trusted_kmers)
            corrected_read[i:i+k] = corrected_kmer
    return ''.join(corrected_read)

# Example trusted k-mers and erroneous read
trusted_kmers = {"ACGCT", "CGTCA", "GTCAA"}
erroneous_read = "ACGCTCGTCC"

# Correct the read
corrected_read = correct_read_using_neighbors(erroneous_read, 5, trusted_kmers)
print("Corrected Read:", corrected_read)
```

In this example, the `correct_read_using_neighbors` function looks for erroneous k-mers in the read and attempts to correct them based on the closest trusted k-mer. It applies the Hamming distance method to minimize differences and ensure that the corrected read aligns with trusted k-mers.

### **Summary**

- **Hamming Distance** is used to find the minimum number of nucleotide changes required to convert an erroneous k-mer into a trusted one.
- **Neighboring K-mers** help detect and correct errors by utilizing overlaps between consecutive k-mers. The correction process ensures the consistency of k-mers within the read and reduces sequencing errors effectively.

## 6. **Correct the Entire Read**

- After correcting the erroneous k-mers, reconstruct the full read. Once the errors have been corrected, you can write the corrected sequence back to a FASTQ or FASTA file.

## 7. **Write Corrected Reads**

- Use Biopython to output the corrected reads to a new file in the FASTA or FASTQ format:

   ```python
   from Bio.SeqRecord import SeqRecord
   from Bio.Seq import Seq
   
   corrected_records = []  # Store corrected reads here

   for record in corrected_reads:
       corrected_seq = Seq(corrected_sequence)  # Corrected sequence string
       new_record = SeqRecord(corrected_seq, id=record.id, description="")
       corrected_records.append(new_record)

   with open("corrected_reads.fastq", "w") as output_handle:
       Bio.SeqIO.write(corrected_records, output_handle, "fastq")
   ```

## 8. **Validation**

- After correction, you might want to validate the reads by mapping them back to a reference genome to ensure the errors have been properly corrected and there are no new introduced errors.

---

## Example Workflow Summary

1. **Load reads** using Biopython.
2. **Build K-mer table** to count occurrences.
3. **Set a threshold** to identify trusted and erroneous k-mers.
4. **Detect and correct errors** based on trusted k-mers.
5. **Output corrected reads** back to FASTA/FASTQ format.

This process ensures that low-frequency errors in Illumina reads are reduced, thereby improving the quality of downstream analyses.
