#!/usr/bin/env python

"""
Identify kmers from Illumina Sequences data to correct errors within the dataset.
Plot the distribution of the kmers on a histogram
"""

# Imports
import sys
import numpy as np

from Bio.Seq import Seq
from Bio import SeqIO

# System Arguments (Format: correct_reads.py <reads to correct> <format> <kmer size> <threshold>)
try:
    FILENAME = sys.argv[1]
    FORMAT = sys.argv[2].lower()
    K_LENGTH = int(sys.argv[3])
    THRESHOLD = int(sys.argv[4])
except ValueError:
    sys.exit("Aborting execution: error parsing system arguments")

def getKmerFrequency():
    # Get Frequency of kmers with a hashmap
