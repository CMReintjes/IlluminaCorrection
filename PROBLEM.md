# Problem

Sequencing errors can have a significant impact on downstream bioinformatics analysis, particularly genome assembly and polymorphism identification. It is recommendevd that errors in input DNA sequences are reduced through trimming of low-quality bases or through an error correction step. One of the simplest ways to correct sequencing errors is to rely on kmer frequency to first identify bases that could be erroneous and second, identify the nucleotide the base should be. In this assignment you will build your own Illumina read correction software. When we sequence a genome for assembly or for SNP discovery, enough sequencing is done to reach some level of genome coverage, X. If `X = 30` this means that on average, each base in the genome has been sequenced 30 times. It also means that each kmer created from the genome should be present an average 30 times. If a sequencing read has an error in it, then any kmers that include that erroneous base will have very low counts (usually 1 or 2) compared to kmers from error-free reads. To correct the base it is simply a matter of changing that base to the other nucleotides and determining if any result in kmers with higher counts.

A strategy to correct reads would look something like this;

1. Count all the kmers in the file of reads you want to correct. This creates the data you’ll need to identify erroneous bases and to correct the base
2. Process each read you want to correct one at a time
3. For each base in a read, check to see if any kmers that include the current base have counts that meet some minimum threshold. Assume that if any kmers with high counts can be generated using the current base then that base is NOT an error.
4. Once you’ve determined if a base in a read is an error, you need to determine which nucleotide that base should be. This can be done by substituting the erroneous base with the other nucleotides and seeing if any of the changes results in high count kmers. If so, then use that base to correct the read.
5. Print the corrected read to a file

## Deliverables

For this assignment you will write your own program to correct Illumina sequencing errors. Your program should do the following;

1. From a file of sequences provided by the user, identify all the kmers in the file and count their occurrence
2. Using the kmer frequencies from step 1, correct the sequencing errors in the file of sequences
3. Plot kmer frequency distribution for the uncorrected and corrected reads and write these to files in ‘PNG’ format
4. Write the corrected sequences to a file in the same format as the original file
5. Your program should use BioPython to parse the sequence file

### Structure

The program will have the following command-line structure

```bash
correct_reads.py <reads to correct> <format> <kmer size> <threshold>
```

`reads to correct` – A file of Illumina reads that are to be corrected.\
`format` – The sequence format of the reads to correct. Only `FASTQ` and `FASTA` are valid entries\
`kmer size` – The kmer size to use for the analysis\
`threshold` – The minimum count a kmer must have to not be considered an error.

#### Example

```bash
correct_reads.py my_file.fastq FASTQ 21 2
```

## What to hand-in

1. The python program `correct_read.py`
2. A k-mer frequence histogram of the `reads_to_correct.fq` data before error correction
3. A k-mer frequence histogram of the `reads_to_correct.fq` data after error correction

## Marking

The program will be marked on the following criteria;

- Completeness of requirements
- Clarity of the Python code (I encourage the use of comments to clarify complex code blocks)
- Does it run using the provided test data
- Does it run using untested data

HINT: Create a very small ‘test’ data set to work with. For example, take a fastq format sequence and
pasted it to a file 5 times and then change a few bases in one of the sequences to simulate errors. Once
your program can correct those errors test it on a larger one.
