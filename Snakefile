import os
import json
import pandas as pd
import csv
import gzip
import yaml


fw_files = glob_wildcards(config["directory"] + "/{prefix}_R1_{suffix}.gz")
files_single = glob_wildcards(config["directory"] + "/{basename}.gz")

if config["paired"]:
	include:
		"rules/paired_end.smk"
elif not config["paired"]:
	include:
		"rules/single_end.smk"

rule all:
	input:
		os.path.join(config["output"], "logs/config_file.yaml"),
		os.path.join(config["output"], "12_report/multiqc_report.html"),
		os.path.join(config["output"], "11_merged/community_and_tax_merged.txt"),
		os.path.join(config["output"], "10_taxonomy/krona_plot.html"),


rule save_config:
	""" Rule to save the config file in the logs directory """
	output:
		config_file_out = os.path.join(config["output"], "logs/config_file.yaml")
	message:
		"Saving config file to log directory"
	run:
		with open(output.config_file_out, "w") as file:
			yaml.dump(config, file)
