#!/usr/bin/env python3

import argparse
import biom
import pandas as pd
from skbio.tree import TreeNode
from skbio.diversity import alpha_diversity
import sys


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Phylogenetic tree in Newick format")
    parser.add_argument(
        "--biom",
        help="Counts per ASV/OTU in BIOM format; names must correspond to leaves in the phylogeny passed to --tree",
    )
    parser.add_argument(
        "--out", help="Path to output file, will be written in TSV format"
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    try:
        args = {
            "tree": snakemake.input["tree"],
            "biom": snakemake.input["biom"],
            "out": snakemake.output[0],
        }
        sys.stderr = open(snakemake.log[0], "w")
    except NameError:
        args = vars(args)
    with open(args["biom"]) as fh:
        tab = biom.parse_table(fh).to_dataframe().transpose()
    tree = TreeNode.read(args["tree"]).root_at_midpoint()
    df = pd.DataFrame(
        {
            "names": tab.index,
            "pd": alpha_diversity("faith_pd", tab, otu_ids=tab.columns, tree=tree),
        }
    ).set_index("names")
    df.to_csv(args["out"], sep="\t")
