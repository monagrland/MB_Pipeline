#!/usr/bin/env python3

import pandas as pd
import sys

if snakemake:
    sys.stderr = open(snakemake.log[0], "w")
    community_table = pd.read_csv(
        snakemake.input["community_table"], sep="\t", index_col=0
    )
    tax_table = pd.read_csv(snakemake.input["tax_table"], sep="\t", index_col=1).iloc[
        :, 1:
    ]
    ASV_with_size_lst = tax_table.index.tolist()
    ASV_sans_size_lst = []
    for ASV in ASV_with_size_lst:
        ASV_sans_size_lst.append(ASV.split(";")[0])
    tax_table.index = ASV_sans_size_lst
    merged_table = pd.merge(
        community_table, tax_table, left_index=True, right_index=True
    )
    merged_table.to_csv(snakemake.output[0], sep="\t")
