#!/usr/bin/env python

"""
Generate test data for Illumina Correction software into specified format
"""

# Import
import os
import sys
import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# System Arguments
try:
    FILENAME = sys.argv[1].lower()
    FORMAT = sys.argv[2].lower()
    SEQNUM = int(sys.argv[3])
    LENGTH = int(sys.argv[4])
    ERRORNUM = int(sys.argv[5])
    try:
        USESEED = sys.argv[6].lower()
    except:
        USESEED = 'n'
except ValueError:
    sys.exit("Aborting execution: error parsing system arguments")



def useSeed(seed):
    # Use seed for random generation
    np.random.seed(seed)


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
    if USESEED == 'y':
        useSeed(seed=42)

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
        description="Generated test sample with "+str(ERRORNUM)+" errors for [illumina correction software]"
    )
    seqObject.append(errorSeq)
    writeSeqFile(seqObject, FILENAME, FORMAT)


main()
