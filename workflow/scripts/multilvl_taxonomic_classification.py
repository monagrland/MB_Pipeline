#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 13:15:15 2021

@author: krueger_l
"""
import subprocess
import pandas as pd
from Bio import SeqIO
import os
import argparse
import re
import csv
import sys
import numpy as np

TAXONOMIC_RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def arg():
    parser = argparse.ArgumentParser("Multilevel direct taxonomic classification")
    parser.add_argument(
        "-d",
        "--db_list",
        nargs="+",
        default=[],
        help="List of paths to databases for the direct taxonomic classification.",
    )
    parser.add_argument("-z", "--ASVs", help="Path to ASVs")
    parser.add_argument(
        "-t",
        "--threshold",
        help="Threshold for the direct classification with usearch_global",
    )
    parser.add_argument(
        "-s",
        "--sintax_cutoff",
        help="Sintax Cutoff for the hierarchical classification",
    )
    parser.add_argument("-o", "--output", help="Prefered name of the output file")
    parser.add_argument("--krona", help="Path to write table to draw krona plot")
    parser.add_argument("--stats_mqc", help="Path to write summary stats for MultiQC")
    parser.add_argument(
        "-n", "--threads", default=1, help="Number of threads to be used"
    )
    parser.add_argument(
        "-p",
        "--hierarchical_db",
        help="Path for the database for the hierarchical classification",
    )
    parser.add_argument(
        "-k",
        "--keep_results",
        help="Keep Intermediate results of the different levels of taxonomic classification or not. Accepts 'True' and 'False'",
        type=bool,
        default=False,
    )
    parser.add_argument("-l", "--log", help="Path to Logfile", required=False)
    return parser.parse_args()


def direct_classification(ASVs, db, threshold, output, threads, run_nr, log):
    """
    Function to start the direct Classification.

    Parameters
    ----------
    ASVs : TYPE
        DESCRIPTION.
    db : str
        path to the database
    threshold : float
        threshold for the direct classification, e.g. 0.97
    output : str
        path to the output file
    threads : int
        number of available threads
    run_nr : int
        run number, multilevel classification allows multiple direct
        classification runs, this counter increases on every level
    log : str
        path to logfile

    Returns
    -------
    None.

    """
    print(f"##### Direct Vsearch Classification level {run_nr} #####")
    output = f"{os.path.splitext(output)[0]}.{run_nr}.direct.txt"
    output_sam = f"{os.path.splitext(output)[0]}.sam"
    cmd = f"vsearch --usearch_global {ASVs} -db {db} --id {threshold} --uc {output} --threads {threads} -samout {output_sam} -maxaccepts 100 --strand both 2>> {log}"
    subprocess.run(cmd, shell=True)
    print(f"##### Finished Direct Vsearch Classification level {run_nr} #####")


def format_dir_classification(ASVs_path, output, run_nr):
    """
    Function to format the output of the direct classification.

    Parameters
    ----------
    ASVs_path : str
        path to the community table
    output : str
        path to the output file
    run_nr : int
        run number, multilevel classification allows multiple direct
        classification runs, this counter increases on every level

    Returns
    -------
    nohit_fasta_path : str
        Path to the fasta file containing all not yet classified sequences,
        used for further classifications
    hit_ASVs_df : pandas DataFrame
        dataframe containing the classified ASVs and the corresponding taxonomy

    """
    tax_file = pd.read_csv(
        f"{os.path.splitext(output)[0]}.{run_nr}.direct.txt", sep="\t", header=None
    )
    nohit_ASV_lst = tax_file.loc[tax_file[9] == "*", 8].tolist()
    ASVs_dict = SeqIO.to_dict(SeqIO.parse(open(ASVs_path), "fasta"))
    fasta_lst = []
    for ASV in nohit_ASV_lst:
        fasta_lst.append(ASVs_dict[ASV])
    nohit_fasta_path = f"{os.path.splitext(output)[0]}_nohit_ASVs.fa"
    SeqIO.write(fasta_lst, nohit_fasta_path, "fasta")
    hit_ASVs_df = tax_file.loc[tax_file[9] != "*"]
    return (nohit_fasta_path, hit_ASVs_df)


def hierarchical_classification(nohit_ASVs, db, output, threads, log):
    """
    Function to start the hierarchical classification

    Parameters
    ----------
    nohit_ASVs : str
        path to ASVs that could not be classified during the multilevel
        direct classification process
    db : str
        path to the database
    output : str
        path to the output file
    threads : int
        number of available threads
    log : str
        path to logfile

    Returns
    -------
    None.

    """
    print("##### Hierarchical Vsearch Classification #####")
    output = os.path.splitext(output)[0] + ".hierarchical.txt"
    cmd = f"vsearch --sintax {nohit_ASVs} -db {db} -tabbedout {output} -threads {threads} -strand plus 2>> {log}"
    subprocess.run(cmd, shell=True)
    print("##### Finished Hierarchical Vsearch Classification #####")


def format_hierarchical_classification(output, threshold):
    """
    Function to format the output of the hierarchical classification

    Parameters
    ----------
    output : str
        path to the output file
    threshold : float
        threshold for the hierarchical classification, e.g. 0.97

    Returns
    -------
    names_df : pandas DataFrame


    """

    # The following part is needed to avoid an error that appears if there are
    # empty cells in a csv.

    path = os.path.splitext(output)[0] + ".hierarchical.txt"
    d = []
    with open(path) as csv_file:
        areader = csv.reader(csv_file)
        max_elems = 0
        for row in areader:
            if max_elems < len(row):
                max_elems = len(row)
        csv_file.seek(0)
        for i, row in enumerate(areader):
            d.append(row + ["" for x in range(max_elems - len(row))])
    tax_df = pd.DataFrame(d)

    first_col = tax_df.iloc[:, 0]

    otu_name_lst = []
    domain_lst = []
    for entry in first_col.tolist():
        split = entry.split()
        otu_name_lst.append(split[0])
        if len(split) > 1:
            domain_lst.append(split[1])
        else:
            domain_lst.append("")

    tax_df.index = otu_name_lst
    tax_df[0] = domain_lst

    last_row = tax_df.iloc[:, 6]

    species_lst = []
    for j in last_row.tolist():
        if len(j) >= 1:
            species_lst.append(j.split("+")[0])
        else:
            species_lst.append("")
    tax_df[6] = species_lst
    tax_df = tax_df.iloc[:, :7]
    names_lst = []
    values_lst = []
    for i, r in tax_df.iterrows():
        names = []
        values = []
        for entry in r.tolist():
            if len(entry) > 0:
                split = entry.rsplit("(", 1)
                names.append(split[0])
                values.append(
                    split[1].replace(")", "").replace("\t", "").replace("+", "")
                )
            else:
                split = ""
                names.append(split)
                values.append(0.00)

        names_lst.append(names)
        values_lst.append(values)
    for i, row in enumerate(values_lst):
        for j, value in enumerate(row):
            if float(value) <= float(threshold):
                names_lst[i][j] = ""
            else:
                names_lst[i][j] == ""

    names_df = pd.DataFrame(names_lst)
    names_df.columns = TAXONOMIC_RANKS
    names_df.insert(0, "ASV", otu_name_lst)
    names_df["Kingdom"] = names_df.loc[:, "Kingdom"].str.replace("d:", "k:")
    return names_df


def get_tax_from_samfile(samfile_path):
    """
    Function to read the samfiles created during the classification,
    checks whether there was more than one hit during direct classification
    and if yes, checks for the lowest common ancestor of all direct hits above
    the specified threshold

    Parameters
    ----------
    samfile_path : str
        path to the samfile

    Returns
    -------
    cleaned_taxonomy_table : pandas DataFrame
        Taxonomy table with LCA instead of direct hits if more than one direct
        hit above the threshold

    """
    if os.path.getsize(samfile_path) == 0:
        print(f"No hits in {samfile_path}.")
        return pd.DataFrame(columns=["ASV"] + TAXONOMIC_RANKS)

    sam_df = pd.read_csv(samfile_path, sep="\t", header=None)
    sam_df = sam_df.iloc[:, [0, 2]]
    # create seperate columns out of the different taxonomic ranks
    tax_lst = sam_df.loc[:, 2].tolist()
    edited_tax_lst = []
    for tax in tax_lst:
        tax_string = re.split(".;tax=", tax)[
            1
        ]  # removes the identifier at the start of the taxonomy values
        tax_string = tax_string.removesuffix(";")
        edited_tax_lst.append(tax_string)

    sam_df.drop(2, axis=1, inplace=True)  # remove taxonomy from table
    sam_df = sam_df.rename({0: "ASV"}, axis="columns")
    tax_lst_of_lst = [tax_str.split(",") for tax_str in edited_tax_lst]
    # The following loop is needed to account for missing taxonomic
    tax_lst_of_lst_with_NA = []
    for tax_lst in tax_lst_of_lst:
        ranks_let = ["k", "p", "c", "o", "f", "g", "s"]
        # taxon elements are formatted like: k:Metazoa
        taxdict = {
            taxon.split(":")[0]: ":".join(taxon.split(":")[1:]) for taxon in tax_lst
        }
        taxonomy_with_NA = [
            (":".join([r, taxdict[r]]) if r in taxdict else r + ":NA")
            for r in ranks_let
        ]
        tax_lst_of_lst_with_NA.append(taxonomy_with_NA)

    sam_df[TAXONOMIC_RANKS] = tax_lst_of_lst_with_NA  # add taxonomy as seperate columns

    unique_ASV_lst = list(
        set(sam_df["ASV"].tolist())
    )  # gets the name of all ASVs and removes duplicates
    taxonomy_lst_of_lst = []
    for ASV in unique_ASV_lst:
        ASV_df = sam_df[
            sam_df.loc[:, "ASV"] == ASV
        ]  # create new dataframe for every ASV
        if len(ASV_df) == 1:
            tax_row = ASV_df.iloc[0].tolist()
            taxonomy_lst_of_lst.append(tax_row)
        else:
            tax_row = [ASV]
            for tax_rank in TAXONOMIC_RANKS:
                unique_tax_set = set(
                    ASV_df.loc[:, tax_rank]
                )  # removes duplicates from column
                if len(unique_tax_set) == 1:
                    tax_row.append(list(unique_tax_set)[0])
                else:
                    tax_row.append("")
            taxonomy_lst_of_lst.append(tax_row)
    cleaned_taxonomy_table = pd.DataFrame(
        taxonomy_lst_of_lst, columns=(["ASV"] + TAXONOMIC_RANKS)
    )
    return cleaned_taxonomy_table


def statistics(log, stats_table_path):
    """
    A Function to create a statistics table containing information about how
    many ASVs could be classified in each classification step.

    Parameters
    ----------
    log : str
        path to logfile
    stats_table_path : str
        path to the output file containing a statistics table

    Returns
    -------
    None.

    """
    with open(log) as f:
        logfile = f.readlines()
    direct_matches_lines = []
    database_lst = []

    for line in logfile:
        if line.startswith("Matching"):
            direct_matches_lines.append(line)
        elif line.startswith("Classified"):
            hierarchical_matches_line = line
        elif line.startswith("Reading"):
            database_lst.append(os.path.basename(line.split()[2]))

    number_of_hits_lst = []
    for i, line_str in enumerate(direct_matches_lines):
        integers_from_line = [int(x) for x in line_str.split() if x.isnumeric()]
        number_of_hits_lst.append(integers_from_line[0])
        if i == 0:
            number_of_sequences = integers_from_line[1]
    integers_from_line_hierarchical = [
        int(x) for x in hierarchical_matches_line.split() if x.isnumeric()
    ]

    classification_type_lst = []
    for i in range(len(number_of_hits_lst)):
        classification_type_lst.append(f"Direct #{i+1}")
    classification_type_lst.append("Hierarchical")
    number_of_hits_lst.append(integers_from_line_hierarchical[0])
    overall_percentage_lst = [
        f"{round((num_hits/number_of_sequences)*100, 2)}%"
        for num_hits in number_of_hits_lst
    ]
    sequences_left_lst = []
    num_seq = number_of_sequences
    for num_hits in number_of_hits_lst:
        sequences_left_lst.append(num_seq)
        num_seq = num_seq - num_hits
    relative_percentage_lst = []
    for j in range(len(number_of_hits_lst)):
        relative_percentage_lst.append(
            f"{round((number_of_hits_lst[j]/sequences_left_lst[j])*100, 2)}%"
        )

    stat_dict = {
        "Classification Type": classification_type_lst,
        "Database": database_lst,
        "Number of hits": number_of_hits_lst,
        "Out of": sequences_left_lst,
        "Relative Percentage": relative_percentage_lst,
        "Overall Percentage": overall_percentage_lst,
    }
    stat_df = pd.DataFrame(stat_dict)
    stat_df.index.name = "Run Nr."
    stat_df.to_csv(stats_table_path, sep="\t")


def main(args):
    output = args["output"]
    logfile = (
        args["log"] if args["log"] is not None else f"{os.path.splitext(output)[0]}.log"
    )
    count = 0
    for db in args["db_list"]:
        count += 1
        direct_classification(
            args["ASVs"],
            db,
            args["threshold"],
            output,
            str(args["threads"]),
            str(count),
            logfile,
        )
        samfile_path = f"{os.path.splitext(output)[0]}.{count}.direct.sam"
        tax_from_samfile_df = get_tax_from_samfile(samfile_path)
        args["ASVs"], hit_ASVs = format_dir_classification(
            args["ASVs"], output, str(count)
        )
        if count == 1:
            taxonomy_sans_multihits = tax_from_samfile_df
            hits_df = hit_ASVs
        else:
            taxonomy_sans_multihits = pd.concat(
                [taxonomy_sans_multihits, tax_from_samfile_df]
            )
            hits_df = pd.concat([hits_df, hit_ASVs])
        hits_df.to_csv(
            os.path.splitext(output)[0] + ".direct.full.txt",
            sep="\t",
            header=None,
            index=None,
        )
    hierarchical_classification(
        args["ASVs"], args["hierarchical_db"], output, args["threads"], logfile
    )
    formatted_hierarchical_classification = format_hierarchical_classification(
        output, args["sintax_cutoff"]
    )
    if len(hits_df) > 0:
        formatted_direct_classification = taxonomy_sans_multihits
        taxonomy_df = pd.concat(
            [formatted_direct_classification, formatted_hierarchical_classification],
            ignore_index=True,
        )
    else:
        taxonomy_df = formatted_hierarchical_classification

    # At the moment, the taxonomy_df dataframe contains cells with the string
    # "", but some scripts cant work with that, thats why we replace it in the
    # next few lines

    replace_values = ["k:", "p:", "c:", "o:", "f:", "g:", "s:"]
    for column, value in zip(TAXONOMIC_RANKS, replace_values):
        taxonomy_df[column] = taxonomy_df[column].replace("", np.nan).fillna(value)

    taxonomy_df.to_csv(output, sep="\t")
    tax_df_for_krona = taxonomy_df.drop(["ASV"], axis=1)
    tax_df_for_krona.to_csv(args["krona"], sep="\t", header=False, index=False)

    statistics(logfile, args["stats_mqc"])

    if not args["keep_results"]:
        out_directory = os.path.dirname(output)
        subprocess.run(f"rm {out_directory}/*.direct.*", shell=True)
        subprocess.run(f"rm {out_directory}/*.hierarchical.*", shell=True)
        subprocess.run(f"rm {out_directory}/*_nohit_*", shell=True)


if __name__ == "__main__":
    try:  # Called from within Snakemake pipeline
        args = {
            "db_list": snakemake.input["direct_dbs"],
            "hierarchical_db": snakemake.input["hierarchical_db"],
            "threshold": snakemake.config["classification_threshold"],
            "keep_results": snakemake.params["keep_results"],
            "sintax_cutoff": snakemake.params["hierarchical_threshold"],
            "ASVs": snakemake.input["ASVs"],
            "output": snakemake.output["base"],
            "krona": snakemake.output["krona"],
            "stats_mqc": snakemake.output["stats_mqc"],
            "threads": snakemake.threads,
            "log": snakemake.log[0],
        }
    except NameError:  # Called from commandline
        args = vars(arg())  # convert to dict
    main(args=args)
