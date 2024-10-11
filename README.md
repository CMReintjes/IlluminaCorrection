# IlluminaCorrection

## Author

Connor Reintjes - 400312785

## Association

BIOTECH 4BI3 Bioinformatics: Assignment 1\
McMaster Univeristy - W Booth School of Engineering Practice and Technology

## Description

Basic kmer based error correction software for illumina file processing. Written for use with command line arguments for running through command line interface. Can take both Fastq and Fasta formats as inputs, and will output the corrected file with the same format.

### Running the script

```bash
correct_reads.py [illumina_file] [format] [kmer_size] [threshold] [optional-flags]
```

#### Mandatory Arguments

`illumina_file` – A file of Illumina reads that are to be corrected.\
`format` – The sequence format of the reads to correct. Only `FASTQ` and `FASTA` are valid entries\
`kmer_size` – The kmer size to use for the analysis\
`threshold` – The minimum count a kmer must have to not be considered an error.

#### Optional Flags

`-v`, `--verbose`: Enables verbose output. Verbose output shows all records being processed, the erroneous kmers and corrected.\
`-s`, `--save-plot`: Autosaves plot to the same directory as the software in .png format using the plot titles as filenames.

## Package Requirements

- Matplotlib
- Biopython
