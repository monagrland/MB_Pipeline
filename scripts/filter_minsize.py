#!/usr/bin/env python3

# Filter Fasta file by minimum size parameter in the header
# Example header:
# > sequence_id;size=1234;
# Optional trailing ;

# usage: python filter_minsize.py input.fasta 8 > output.fasta

import re
import sys
from Bio import SeqIO
from sys import stdout

try: # Run within Snakemake pipeline
    infile = snakemake.input[0]
    minsize = int(snakemake.wildcards["minsize"])
    outfile = snakemake.output[0]
except NameError: # Run from command line
    infile = sys.argv[1]
    minsize = int(sys.argv[2])
    outfile = sys.argv[3]

with open(outfile, "w") as fh:
    for seq in SeqIO.parse(infile, "fasta"):
        size = int(re.search(r"size=(\d+)", seq.id).groups()[0])
        if size >= minsize:
            SeqIO.write(seq, fh, "fasta-2line") # no line wrapping of sequence
