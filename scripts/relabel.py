#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 10:37:07 2023

@author: krueger_l
"""
import os
import gzip
import argparse


def arg():
    parser = argparse.ArgumentParser("Relabel FASTQ Headers")
    parser.add_argument(
        "-i",
        "--input",
        help="Path to the input file (has to be .fastq.gz)",
    )
    parser.add_argument(
        "-o", "--output", help="Path to the output file (has to be .fastq)"
    )
    return parser.parse_args()


def update_headers(input_file, output_file):
    # Extract the filename without the extension
    basename = os.path.basename(
        os.path.splitext(os.path.splitext(input_file)[0].rstrip("/"))[0]
    )

    # Open the input gzipped fastq file for reading
    with gzip.open(input_file, "rt") as infile:
        # Open the output gzipped fastq file for writing
        with open(output_file, "wt") as outfile:
            line_count = 0
            header_count = 0

            for line in infile:
                line_count += 1

                # Identify the header lines every 4th line starting from line 1
                if line_count % 4 == 1:
                    header_count += 1
                    new_header = f"@{basename}_{header_count};sample={basename}"
                    outfile.write(f"{new_header}\n")
                else:
                    # Write the rest of the lines as they are
                    outfile.write(line)


def main():
    args = arg()
    input_file = args.input
    output_file = args.output
    update_headers(input_file, output_file)


if __name__ == "__main__":
    main()
