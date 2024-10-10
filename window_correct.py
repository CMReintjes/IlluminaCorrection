from tracemalloc import start
from Bio.Seq import Seq

def fix_window(window, k, threshold, frequency, reverse=False):
    "Fix the sliding window"
    corrected = list(window)
    all_kmers = True
    if reverse:
        pos = len(window)-k
        start_pos = len(window)-k
        end_pos = -1
        step = -1
        original = ''.join(window[start_pos:start_pos+k+1])
    else:
        pos = k-1
        start_pos = 0
        end_pos = len(window)-k+1
        step = 1
        original = ''.join(window[start_pos:start_pos+k])
    max_freq = frequency[original]
    bases = list("ATGC")
    for b in bases:
        all_kmers = True
        window[pos] = b
        for n in range(start_pos, end_pos, step):
            kmer = "".join(window[n:n+k])
            if kmer not in frequency:
                all_kmers = False
                break

        kmer = list(original)
        kmer[pos] = b
        kmer = ''.join(kmer)
        if all_kmers and kmer in frequency:
            if frequency[kmer] > max_freq and frequency[kmer] > threshold:
                corrected[pos] = b

    return list(corrected[start_pos:start_pos+k])
