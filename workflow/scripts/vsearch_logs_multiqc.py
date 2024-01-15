#!/usr/bin/env python3

# Parse VSEARCH log messages for MultiQC

# import argparse
import re
from os.path import basename, splitext
import json
import argparse

LOGFILE_REGEX = {
    "fastq_mergepairs": {
        "vsearch_version": r"^vsearch ([^,]+),",
        "pairs": r"^\s+(\d+)\s+Pairs$",
        "merged": r"^\s+(\d+)\s+Merged",
        "notmerged": r"^\s+(\d+)\s+Not merged",
    },
    "fastq_filter": {
        "vsearch_version": r"^vsearch ([^,]+),",
        "kept": r"(\d+) sequences kept \(of which \d+ truncated\), \d+ sequences discarded.",
        "kept_truncated": r"\d+ sequences kept \(of which (\d+) truncated\), \d+ sequences discarded.",
        "discarded": r"\d+ sequences kept \(of which \d+ truncated\), (\d+) sequences discarded.",
    },
    "derep_fulllength": {
        "vsearch_version": r"^vsearch ([^,]+),",
        "sequences": r"\d+ nt in (\d+) seqs,",
        "unique": r"(\d+) unique sequences",
    },
}


def onefile(infile, regex):
    """Parse single log file from VSEARCH"""
    out = {}
    with open(infile, "r") as fh:
        for line in fh:
            s = {key: re.match(regex[key], line.rstrip()) for key in regex}
            for key in s:
                if s[key]:
                    out[key] = s[key].groups()[0]
    return out


def process_record(rec, which):
    if which == "fastq_mergepairs":
        for key in ["pairs", "merged", "notmerged"]:
            rec[key] = int(rec[key])
        if rec["merged"] + rec["notmerged"] != rec["pairs"]:
            print("Warning: numbers of reads do not match")
        out_keys = ["merged", "notmerged"]
    elif which == "fastq_filter":
        for key in ["kept", "kept_truncated", "discarded"]:
            rec[key] = int(rec[key])
        rec["kept_ok"] = rec["kept"] - rec["kept_truncated"]
        out_keys = ["kept_ok", "kept_truncated", "discarded"]
    elif which == "derep_fulllength":
        for key in ["sequences", "unique"]:
            rec[key] = int(rec[key])
        rec["replicate"] = rec["sequences"] - rec["unique"]
        out_keys = ["unique", "replicate"]
    return {key: rec[key] for key in out_keys}


def parsefiles(filelist, which, sample_from_filename=True):
    out = {}
    for logfile in filelist:
        if sample_from_filename:
            sample = splitext(basename(logfile))[0]
        else:
            sample = logfile
        out[sample] = process_record(onefile(logfile, LOGFILE_REGEX[which]), which)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--format", type=str, help="VSEARCH subcommand that produced the file"
    )
    parser.add_argument("--files", nargs="+")
    args = parser.parse_args()

    data = parsefiles(args.files, args.format)
    if args.format == "fastq_mergepairs":
        out = {
            "id": "vsearch_fastq_mergepairs",
            "section_name": "Paired end read merging",
            "plot_type": "bargraph",
            "data": data,
        }
    elif args.format == "fastq_filter":
        out = {
            "id": "vsearch_fastq_filter",
            "section_name": "Quality filtering of merged/single reads",
            "plot_type": "bargraph",
            "data": data,
        }
    elif args.format == "derep_fulllength":
        out = {
            "id": "vsearch_derep_fulllength",
            "section_name": "Dereplication of individual read files",
            "plot_type": "bargraph",
            "data": data,
        }
    print(json.dumps(out, indent=4))
