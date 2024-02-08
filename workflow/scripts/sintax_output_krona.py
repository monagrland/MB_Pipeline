#!/usr/bin/env python3

import sys
import re

ranks = ["k", "p", "c", "o", "f", "g", "s"]

if __name__ == "__main__":
    try:
        infile = snakemake.input[0]
        outfile_krona = snakemake.output["krona"]
        outfile_tab = snakemake.output["tax_table"]
        cutoff = snakemake.params["cutoff"]
    except NameError:
        infile = sys.argv[1]
        outfile_krona = sys.argv[2]
        outfile_tab = sys.argv[3]
        cutoff = float(sys.argv[4])

    fh_in = open(infile, "r")
    fh_out_krona = open(outfile_krona, "w")
    fh_out_tab = open(outfile_tab, "w")
    fh_out_tab.write(
        "\t".join(
            [
                "",
                "ASV",
                "Kingdom",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
            ]
        )
        + "\n"
    )
    counter = 0
    for line in fh_in:
        spl = line.rstrip().split("\t")
        if len(spl) >= 2:
            tax = spl[1]
            tax = [re.search(r"(.):(.+)\((.+)\)", i).groups() for i in tax.split(",")]
            tax = {i[0]: i[1] for i in tax if float(i[2]) >= cutoff}
            tax = [
                rank + ":" + tax[rank] if rank in tax else rank + ":" for rank in ranks
            ]
            fh_out_krona.write("\t".join(tax) + "\n")
        elif len(spl) == 1:
            tax = [rank + ":" for rank in ranks]
        fh_out_tab.write(
            "\t".join(
                [
                    str(counter),
                    spl[0],
                ]
                + tax
            )
            + "\n"
        )
        counter += 1
    fh_in.close()
    fh_out_krona.close()
    fh_out_tab.close()
