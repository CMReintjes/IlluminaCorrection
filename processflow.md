# Process Flow

## 1. **Calculate k-mers and their Frequency (k-counts)**

- The DNA sequence is split into small overlapping pieces called _k-mers_ of length "k."
- Each _k-mer_ has a count (or frequency), called the k-count, representing how often that _k-mer_ appears in the dataset.
- If a _k-mer_ has a **high k-count**, it’s considered a "solid" (likely correct) _k-mer_. If it has a **low k-count**, it’s considered "weak" (likely to have errors).
- Because sequencing data can be huge, a hashed map (a fast way of storing data) is used to store each _k-mer_ and its positions in the read.

## 2. **Establish the Error Threshold**

- This step separates solid _k-mers_ (correct) from weak _k-mers_ (incorrect).
- The error threshold is a k-count value that distinguishes between solid and weak _k-mers_.
- To find this threshold, the frequency of each k-count is examined. Typically, the first minimum in the frequency distribution of k-counts is used to set the threshold.
- In simpler terms, this means finding the point where the count of _k-mers_ shifts from being common (solid) to rare (weak).
- If there is a long section in the frequency distribution with consecutive zeros (no _k-mers_), this is used to set the error threshold.

## 3. **Identify Error Regions in the Sequence**

- The sequence is scanned for regions where _k-mers_ have k-counts below the error threshold. These regions are flagged as potential error zones.
- The flagged regions are grouped into clusters using a technique called variable bandwidth mean-shift clustering (to identify error areas more accurately).
- After identifying these error sections, they can be extended in both directions to include adjacent _k-mers_ that might belong to the same error cluster.

## 4. **Correct Errors in the Error Sections**

- Once error regions are identified, the system attempts to correct them.
- It analyzes the surrounding solid _k-mers_ to suggest a correction, replacing the weak _k-mers_ with more likely, high-frequency ones.
- This process improves the accuracy of the DNA sequence data by fixing these error-prone regions.

## Summary

- **Phase 1:** Calculate how often each _k-mer_ appears in the sequence.
- **Phase 2:** Find a threshold that separates common (solid) _k-mers_ from rare (weak) ones.
- **Phase 3:** Identify regions in the DNA sequence where there are likely errors based on the threshold.
- **Phase 4:** Correct the errors using information from the solid (correct) _k-mers_ around them.

This method makes use of the natural repetition of DNA sequences to distinguish between real sequence data and sequencing errors.
