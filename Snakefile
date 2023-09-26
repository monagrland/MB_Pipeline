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


rule all:
	input:
		"logs/config_file.yaml",
		"12_report/multiqc_report.html",
		"11_merged/community_and_tax_merged.txt",
		"10_taxonomy/krona_plot.html",
		# expand("07_ASVs/ASVs_{alpha}_entropy_values.csv", alpha=[1,2,3,4,5,6,7,8,9,10]),


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
