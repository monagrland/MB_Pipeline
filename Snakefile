import os
import json
import pandas as pd
import csv
import gzip
import yaml
from snakemake.utils import validate

validate(config, "config/config.schema.yaml")

workdir : config['workdir']
reads_df = pd.read_table(config['reads_table'], sep="\t").set_index("sample", drop=False)
samples = reads_df['sample'].drop_duplicates()

if config["paired"]:
	include:
		"rules/paired_end.smk"
elif not config["paired"]:
	include:
		"rules/single_end.smk"

include: "rules/common.smk"
include: "rules/coding.smk"

if config["protein_coding"]:
	screening="no_pseudogenes"
else:
	screening="no_chimeras"

rule all:
	input:
		"logs/config_file.yaml",
		expand(
			[
				"10_taxonomy/krona_plot.{method}.{screening}.html",
				"11_merged/community_and_tax_merged.{method}.{screening}.txt",
				"12_report/multiqc_report.{method}.{screening}.html",
				"13_phylogeny/ASVs_{method}.{screening}.{treeprog}.faiths_pd.tsv",
			],
			method=config['denoising']['method'],
			screening=screening,
			treeprog=config['treeprog'],
		),

rule save_config:
	""" Rule to save the config file in the logs directory """
	output:
		config_file_out = "logs/config_file.yaml"
	message:
		"Saving config file to log directory"
	run:
		with open(output.config_file_out, "w") as file:
			yaml.dump(config, file)

# vim: set noexpandtab:
