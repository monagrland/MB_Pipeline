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
		# "logs/config_file.yaml",
		# "12_report/multiqc_report.html",
		# "11_merged/community_and_tax_merged.txt",
		# "10_taxonomy/krona_plot.html",
		"diagnostics/entropy_ratio_denoising_plot.png"


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
