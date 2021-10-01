import os
import json
import pandas as pd
import csv
import gzip
import shutil


fw_files = glob_wildcards(config["directory"] + "/{prefix}_R1_{suffix}.gz")


def output_lst():
	"""
	Function to specify the wanted output.
	Especially used during debugging or when you dont want specific output files
	like the MultiQC report or the krona plot
	"""
	out_lst = [config["output"] + "/09_OTU_table/otu_table.txt"]
	if config["mqc_report"]:
		out_lst.append(config["output"] + "/10_report/multiqc_report.html")
		print("Added MultiQC report to Output")
	if config["taxonomy"]:
		out_lst.append(config["output"] + "/11_taxonomy/krona_plot.html")
		# out_lst.append(config["output"] + "/12_merged/otu_and_tax_merged.txt")
		print("Added Taxonomy to Output")
	return(out_lst)

rule all:
	input:
		output_lst()


# rule remove_short_reads:
# 	input:
# 		input_fw = config["directory"] + "/{prefix}_R1_{suffix}.gz",
# 		input_rv = config["directory"] + "/{prefix}_R2_{suffix}.gz"
# 	output:
# 		output_fw = config["output"] + "/00_wo_short_reads/{prefix}_R1_{suffix}.gz",
# 		output_rv = config["output"] + "/00_wo_short_reads/{prefix}_R2_{suffix}.gz"
# 	params:
# 		length_th = int(config["00_length_threshold"])
# 	message:
# 		"Removing files with less than {params} entries or where one of two files is empty"
# 	run:
# 		fw_lst = []
# 		rv_lst = []
# 		with gzip.open(input.input_fw, "rb") as f:
# 		    for line in f:
# 		        fw_lst.append(line)
# 		with gzip.open(input.input_rv, "rb") as g:
# 		    for line in g:
# 		        rv_lst.append(line)
# 		if len(fw_lst) > 100 and len(rv_lst) > 100:
# 			if len(fw_lst) > 0 and len(rv_lst) > 0:
# 				shutil.copyfile(input.input_fw, output.output_fw)
# 				shutil.copyfile(input.input_rv, output.output_rv)
# 			else:
# 				print("Skipping " + input.input_fw + " and " + input.input_rv + "because at least one of the files is empty")
# 		else:
# 			print("Skipping " + input.input_fw + " and " + input.input_rv + "because at least one of the files has less than" + str(params.length_th) + "entries")



rule cutadapt:
	""" Rule to remove the Adapter Sequences from the reads """
	input:
		input_fw = config["directory"] + "/{prefix}_R1_{suffix}.gz",
		input_rv = config["directory"] + "/{prefix}_R2_{suffix}.gz"
	output:
		output_fw = config["output"] + "/01_trimmed_data/{prefix}_R1_{suffix}.gz",
		output_rv = config["output"] + "/01_trimmed_data/{prefix}_R2_{suffix}.gz"
	params:
		options = " ".join(config["01_adapter_trimming_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
	conda:
		os.path.join(workflow.basedir, "envs/mb_cutadapt.yaml")
	message:
		"Executing adaptertrimming for {params.filename_fw} and {params.filename_rv}"
	log:
		config["output"] + "/logs/01_cutadapt/{prefix}_{suffix}.txt"
	shell:
		"cutadapt {params.options} -o {output.output_fw} -p {output.output_rv} {input.input_fw} {input.input_rv} &>>  {log}"

rule trimmomatic:
	"""rule to filter reads"""
	input:
		input_fw = config["output"] + "/01_trimmed_data/{prefix}_R1_{suffix}.gz",
		input_rv = config["output"] + "/01_trimmed_data/{prefix}_R2_{suffix}.gz"
	output:
		output_fw = config["output"] + "/02_filtered_data/{prefix}_R1_{suffix}.gz",
		output_rv = config["output"] + "/02_filtered_data/{prefix}_R2_{suffix}.gz"
	params:
		options = " ".join(config["02_trimmomatic_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
	conda:
		os.path.join(workflow.basedir, "envs/mb_trimmomatic.yaml")
	message:
		"Removing low quality reads for {params.filename_fw} and {params.filename_rv}"
	log:
		config["output"] + "/logs/02_trimmomatic/{prefix}_{suffix}.txt"
	shell:
		"trimmomatic PE {input.input_fw} {input.input_rv} {output.output_fw} /dev/null {output.output_rv} /dev/null {params.options} &>> {log}"

rule merge:
	""" Rule to merge paired end reads to a single file"""
	input:
		input_fw = config["output"] + "/02_filtered_data/{prefix}_R1_{suffix}.gz",
		input_rv = config["output"] + "/02_filtered_data/{prefix}_R2_{suffix}.gz"
	output:
		config["output"] + "/03_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	params:
		options = " ".join(config["03_vsearch_merge_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
		basename = "{prefix}"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Merging paired end reads for {params.filename_fw} and {params.filename_rv}"
	log:
		config["output"] + "/logs/03_merging/{prefix}_{suffix}.txt"
	shell:
		"vsearch --fastq_mergepairs {input.input_fw} --reverse {input.input_rv} --fastaout {output} {params.options} --relabel {params.basename}_ --label_suffix \;sample={params.basename} &>> {log}"

rule dereplicate:
	""" rule to remove duplicates """
	input:
		config["output"] + "/03_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	output:
		config["output"] + "/04_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	params:
		filename = "{prefix}_{suffix}.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing redundant reads for {params.filename}"
	log:
		config["output"] + "/logs/04_dereplicate/{prefix}_" + os.path.splitext("{suffix}")[0] + ".txt"
	shell:
		"vsearch --derep_fulllength {input} --output {output} --sizeout &>> {log}"

rule concatenate:
	""" rule to concatenate all files into one """
	input:
		expand((config["output"] + "/04_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"), zip, prefix = fw_files.prefix, suffix = fw_files.suffix)
	output:
		config["output"] + "/05_concatenated_data/all_reads.fasta"
	params:
		derep_dir = config["output"] + "/04_derep_data/"
	message:
		"Concatenating all reads"
	log:
		config["output"] + "/logs/05_concatenate/all_reads.txt"
	shell:
		"cat {params.derep_dir}*.fasta > {output}"

rule dereplicate_2:
	""" remove the duplicates inside the single file """
	input:
		config["output"] + "/05_concatenated_data/all_reads.fasta"
	output:
		config["output"] + "/06_derep_data/unique_reads.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing redundant reads for all reads"
	log:
		config["output"] + "/logs/06_dereplicate/all_reads.txt"
	shell:
		"vsearch --derep_fulllength {input} --output {output} --sizein --sizeout &>> {log}"

rule generate_zOTUs:
	input:
		config["output"] + "/06_derep_data/unique_reads.fasta"
	output:
		config["output"] + "/07_zOTUs/zOTUs.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Generating zOTUs"
	log:
		config["output"] + "/logs/07_zOTUs/all_reads.txt"
	shell:
		"vsearch --cluster_unoise {input} --centroids {output} --relabel zOTU --sizein --sizeout &>>{log}"

rule remove_chimeras:
	input:
		config["output"] + "/07_zOTUs/zOTUs.fasta"
	output:
		config["output"] + "/08_zOTUs_nonchimeras/zOTUs_nonchimeras.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing Chimeras"
	log:
		config["output"] + "/logs/08_zOTUs_nonchimeras/all_reads.txt"
	shell:
		"vsearch -uchime3_denovo {input} --nonchimeras {output} --sizein --xsize &>> {log}"

rule generate_OTU_table:
	input:
		search = config["output"] + "/06_derep_data/unique_reads.fasta",
		db = config["output"] + "/08_zOTUs_nonchimeras/zOTUs_nonchimeras.fasta"
	output:
		otu_table = config["output"] + "/09_OTU_table/otu_table.txt",
	params:
		options = " ".join(config["09_otu_table_options"]),
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Generating OTU table"
	log:
		config["output"] + "/logs/09_OTU_table/all_reads.txt"
	shell:
		"vsearch --usearch_global {input.search} --db {input.db} {params.options} --otutabout {output.otu_table} &>> {log}"

rule generate_report:
	input:
		otu_table = config["output"] + "/09_OTU_table/otu_table.txt",
		custom_mqc_config = "multiqc_config.yaml"
	output:
		config["output"] + "/10_report/multiqc_report.html"
	conda:
		os.path.join(workflow.basedir, "envs/mb_multiqc.yaml")
	params:
		output_dir = config["output"] + "/10_report/",
		log_dir = config["output"] + "/logs"
	message:
		"Generating MultiQC report"
	log:
		config["output"] + "/logs/10_MultiQC/multiqc.txt"
	shell:
		"multiqc {params.log_dir} -o {params.output_dir} --config {input.custom_mqc_config}"

rule taxonomy:
	input:
		zOTUs = config["output"] + "/08_zOTUs_nonchimeras/zOTUs_nonchimeras.fasta"
	output:
		base = config["output"] + "/11_taxonomy/taxonomy.txt",
		plot = config["output"] + "/11_taxonomy/taxonomy.krona.txt"
	params:
		direct_db_lst = config["11_direct_dbs"],
		hierarchical_db = config["11_hierarchical_db"],
		threshold = config["11_threshold"],
		keep_results = "False"
	threads:
		config["threads"]
	message:
		"Starting Multilevel Taxonomic Classification"
	conda:
		os.path.join(workflow.basedir, "envs/mb_taxonomy.yaml")
	shell:
		"python3 scripts/multilvl_taxonomic_classification.py -d {params.direct_db_lst} -z {input.zOTUs} -t {params.threshold} -o {output.base} -n {threads} -p {params.hierarchical_db} -k {params.keep_results}"


rule krona:
	input:
		config["output"] + "/11_taxonomy/taxonomy.krona.txt"
	output:
		config["output"] + "/11_taxonomy/krona_plot.html"
	conda:
		os.path.join(workflow.basedir, "envs/mb_krona.yaml")
	message:
		"Creating Krona Plot"
	shell:
		"ktImportText -q {input} -o {output}"
