import os
import json
import pandas as pd
import csv
import gzip
import yaml

workdir : config['output']
fw_files = glob_wildcards(config["directory"] + "/{prefix}_R1_{suffix}.gz")
files_single = glob_wildcards(config["directory"] + "/{basename}.gz")

if config["paired"]:
	include:
		"rules/paired_end.smk"
elif not config["paired"]:
	include:
		"rules/single_end.smk"

include: "rules/common.smk"
include: "rules/common_coding.smk"

if config["protein_coding"]:
	screening="no_pseudogenes"
else:
	screening="no_chimeras"

rule all:
	input: # pseudorule for protein coding sequences
		"logs/config_file.yaml",
		# "12_report/multiqc_report.html",
		# "11_merged/community_and_tax_merged.txt",
		# "10_taxonomy/krona_plot.html",
		# expand("08_ASVs_screened/ASVs_{method}.no_pseudogenes.fasta", method=config['denoising']['method'])
		expand("10_taxonomy/krona_plot.{method}.{screening}.html", method=config['denoising']['method'], screening=screening)

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
