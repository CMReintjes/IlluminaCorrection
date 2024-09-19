#!/usr/bin/env python

"""
Generate test data for Illumina Correction software into specified format
"""

# Import
import os
import sys
import argparse
import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# System Arguments
parser = argparse.ArgumentParser()
parser.add_argument('file_output', type=str,
                    help='Name of output file', required=True)
parser.add_argument('-t', '--type', type=str, default='fasta', choices=[
                    'fasta', 'fastq'], help='Specify file format (Default: FASTA)', required=True)
parser.add_argument('-l', '--length', type=int, default=100,
                    help='Length of nucleotide sequences to generate', required=True)
parser.add_argument('-e', '--errors', type=int, default=1,
                    help='Number of errors to generate', required=True)
parser.add_argument('-n', '--number', type=int, default=5,
                    help='Number of sequences to generate', required=True)
parser.add_argument('-s', '--seed', type=int, default=42,
                    help='Use seed to generate tests')

args = parser.parse_args()
FILENAME = args.file_output
FORMAT = args.type
SEQNUM = args.number
LENGTH = args.length
ERRORNUM = args.errors

if args.seed:
    np.random.seed(args.seed)
    

def genSequence(seqLen):
    # Generate a random nucleotide sequence of length <seqLen> and return
    sequence = np.random.choice(list("ACTG"), seqLen)
    return sequence


def genRandomError(seqList, n):
    # Generate a set number <n> of random positions within the nucleotide sequence to have an error
    errorPos = np.random.choice(len(seqList), n)
    for i in range(0, len(errorPos)):
        # Loop through random nucleotide selection if it matches the nucleotide at current position
        nucDiff = True
        while nucDiff:
            # Generate random nucleotide
            errorNuc = np.random.choice(list("ACTG"))
            if seqList[errorPos[i]] != errorNuc[0]:
                break
        seqList[errorPos[i]] = errorNuc[0]
    return seqList


def writeSeqFile(seqObj, fileName, fileType):
    # Write the sequences out into a file with name <fileName>
    SeqIO.write(seqObj, fileName, fileType)
    print("Successful")


def main():
    global FILENAME, FORMAT, SEQNUM, LENGTH, ERRORNUM, USESEED
    seqObject = []

    # Generate the master sequence and add to SeqObject
    masterSequenceList = genSequence(LENGTH)
    masterSeq = ''.join(masterSequenceList)
    for i in range(1, SEQNUM):
        tempSeq = SeqRecord(
            Seq(masterSeq),
            id="illumina|test|seq"+str(i),
            description="Generated test sample for [illumina correction software]"
        )
        seqObject.append(tempSeq)

    # Use the master sequence to generate error positions and add to SeqObject
    errorSequenceList = genRandomError(masterSequenceList, ERRORNUM)
    errorSequence = ''.join(errorSequenceList)
    errorSeq = SeqRecord(
        Seq(errorSequence),
        id="illumina|test|seq"+str(0),
        description="Generated test sample with " +
        str(ERRORNUM)+" errors for [illumina correction software]"
    )
    seqObject.append(errorSeq)
    writeSeqFile(seqObject, FILENAME, FORMAT)


main()
