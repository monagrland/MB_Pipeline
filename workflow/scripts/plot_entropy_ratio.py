#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import os.path
import sys
import argparse


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", help="List of _entropy_values.csv files produced by DnoisE with the -g option, comma-separated")
    parser.add_argument("--alphas", help="Corresponding alpha values for each input file, comma-separated")
    parser.add_argument("--output", default="test.denoised_entropy_ratio.png", help="Output file path")
    return parser.parse_args()


if __name__ == "__main__":
    try: # called from Snakemake pipeline
        args = {
            'inputs' : snakemake.params['inputs'],
            'alphas' : snakemake.params['alphas'],
            'output' : snakemake.output[0],
        }
    except NameError: # called from command line
        args = vars(args())

    print(args['inputs'])
    print(args['alphas'])

    files = args['inputs'].split(",")
    alphas = [int(i) for i in args['alphas'].split(",")]
    dfs = []
    for fn,alpha in zip(files, alphas):
        df = pd.read_csv(fn, sep=",")
        df['alpha'] = alpha
        dfs.append(df)
    df = pd.concat(dfs)
    df['e2e3_ratio'] = df['e2']/df['e3']

    num_lengths = len(df['seq_length'].value_counts())
    if num_lengths > 1:
        print("Multiple sequence lengths, picking the length representing the most sequences")
    lenmax = list(df.groupby("seq_length")[["total_count"]].sum().idxmax())[0]
    fig, axs = plt.subplots(2,1, sharex=True)
    axs[0].scatter(
        df.query(f"seq_length == {lenmax}")["alpha"],
        df.query(f"seq_length == {lenmax}")["e2e3_ratio"],
    )
    axs[0].set_xlim(max(alphas)+0.5, min(alphas)-0.5)
    axs[0].set_ylabel("Entropy ratio")
    axs[1].scatter(
        df.query(f"seq_length == {lenmax}")["alpha"],
        df.query(f"seq_length == {lenmax}")["total_seqs"],
    )
    axs[1].set_xlim(max(alphas)+0.5, min(alphas)-0.5)
    axs[1].set_ylabel("Denoised sequences")
    fig.suptitle(f"Sequence length {lenmax}")
    fig.tight_layout()
    fig.savefig(args['output'])
