import os
import json
import pandas as pd
import csv
import gzip
import shutil


fw_files = glob_wildcards(config["directory"] + "/{prefix}_R1_{suffix}.gz")

def create_gzip(gz_file_name):
	"""
	Function to create a gzip archive with an empty file. It is necessary to
	create empty archives because Snakemake expects the files to exist,
	even if they are empty.
	"""
	unzipped_name = os.path.splitext(gz_file_name)[0]
	os.system("touch " + unzipped_name)
	with open(unzipped_name, "rb") as f_in:
		with gzip.open(gz_file_name, "wb") as f_out:
			shutil.copyfileobj(f_in, f_out)
	os.system("rm " + unzipped_name)

rule all:
	input:
		config["output"] + "/11_report/multiqc_report.html",
		config["output"] + "/12_taxonomy/taxonomy.txt",
		config["output"] + "/12_taxonomy/krona_plot.html",

rule remove_short_files:
	input:
		input_fw = config["directory"] + "/{prefix}_R1_{suffix}.gz",
		input_rv = config["directory"] + "/{prefix}_R2_{suffix}.gz"
	output:
		output_fw = config["output"] + "/01_sans_short_reads/{prefix}_R1_{suffix}.gz",
		output_rv = config["output"] + "/01_sans_short_reads/{prefix}_R2_{suffix}.gz"
	params:
		length_th = int(config["min_file_length_threshold"])
	message:
		"Removing files with less than {params} entries or where one of two files is empty"
	run:
		fw_lst = []
		rv_lst = []
		with gzip.open(input.input_fw, "rb") as f:
		    for line in f:
		        fw_lst.append(line)
		with gzip.open(input.input_rv, "rb") as g:
		    for line in g:
		        rv_lst.append(line)
		if len(fw_lst) > params.length_th or len(rv_lst) > params.length_th:
			if len(fw_lst) > 0 and len(rv_lst) > 0:
				shutil.copyfile(input.input_fw, output.output_fw)
				shutil.copyfile(input.input_rv, output.output_rv)
			else:
				create_gzip(output.output_fw)
				create_gzip(output.output_rv)
				print("Skipping " + os.path.basename(input.input_fw) + " and " + os.path.basename(input.input_rv) + " because at least one of the files is empty")
		else:
			create_gzip(output.output_fw)
			create_gzip(output.output_rv)
			print("Skipping " + os.path.basename(input.input_fw) + " and " + os.path.basename(input.input_rv) + " because at least one of the files has less than " + str(params.length_th) + " entries")



rule cutadapt:
	""" Rule to remove the Adapter Sequences from the reads """
	input:
		input_fw = config["output"] + "/01_sans_short_reads/{prefix}_R1_{suffix}.gz",
		input_rv = config["output"] + "/01_sans_short_reads/{prefix}_R2_{suffix}.gz"
	output:
		output_fw = config["output"] + "/02_trimmed_data/{prefix}_R1_{suffix}.gz",
		output_rv = config["output"] + "/02_trimmed_data/{prefix}_R2_{suffix}.gz"
	params:
		options = " ".join(config["adapter_trimming_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
	conda:
		os.path.join(workflow.basedir, "envs/mb_cutadapt.yaml")
	message:
		"Executing adaptertrimming for {params.filename_fw} and {params.filename_rv}"
	log:
		config["output"] + "/logs/02_cutadapt/{prefix}_{suffix}.txt"
	shell:
		"cutadapt {params.options} -o {output.output_fw} -p {output.output_rv} {input.input_fw} {input.input_rv} &>>  {log}"

rule trimmomatic:
	"""rule to filter reads"""
	input:
		input_fw = config["output"] + "/02_trimmed_data/{prefix}_R1_{suffix}.gz",
		input_rv = config["output"] + "/02_trimmed_data/{prefix}_R2_{suffix}.gz"
	output:
		output_fw = config["output"] + "/03_filtered_data/{prefix}_R1_{suffix}.gz",
		output_rv = config["output"] + "/03_filtered_data/{prefix}_R2_{suffix}.gz"
	params:
		options = " ".join(config["trimmomatic_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
	conda:
		os.path.join(workflow.basedir, "envs/mb_trimmomatic.yaml")
	message:
		"Removing low quality reads for {params.filename_fw} and {params.filename_rv}"
	log:
		config["output"] + "/logs/03_trimmomatic/{prefix}_{suffix}.txt"
	shell:
		"trimmomatic PE {input.input_fw} {input.input_rv} {output.output_fw} /dev/null {output.output_rv} /dev/null {params.options} &>> {log}"

rule merge:
	""" Rule to merge paired end reads to a single file"""
	input:
		input_fw = config["output"] + "/03_filtered_data/{prefix}_R1_{suffix}.gz",
		input_rv = config["output"] + "/03_filtered_data/{prefix}_R2_{suffix}.gz"
	output:
		config["output"] + "/04_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	params:
		options = " ".join(config["vsearch_merge_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
		basename = "{prefix}"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Merging paired end reads for {params.filename_fw} and {params.filename_rv}"
	log:
		config["output"] + "/logs/04_merging/{prefix}_{suffix}.txt"
	shell:
		"vsearch --fastq_mergepairs {input.input_fw} --reverse {input.input_rv} --fastaout {output} {params.options} --relabel {params.basename}_ --label_suffix \;sample={params.basename} &>> {log}"

rule dereplicate:
	""" rule to remove duplicates """
	input:
		config["output"] + "/04_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	output:
		config["output"] + "/05_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	params:
		filename = "{prefix}_{suffix}.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing redundant reads for {params.filename}"
	log:
		config["output"] + "/logs/05_dereplicate/{prefix}_" + os.path.splitext("{suffix}")[0] + ".txt"
	shell:
		"vsearch --derep_fulllength {input} --output {output} --strand plus --sizeout &>> {log}"

rule concatenate:
	""" rule to concatenate all files into one """
	input:
		expand((config["output"] + "/05_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"), zip, prefix = fw_files.prefix, suffix = fw_files.suffix)
	output:
		config["output"] + "/06_concatenated_data/all_reads.fasta"
	params:
		derep_dir = config["output"] + "/05_derep_data/"
	message:
		"Concatenating all reads"
	log:
		config["output"] + "/logs/06_concatenate/all_reads.txt"
	shell:
		"cat {params.derep_dir}*.fasta > {output}"

rule dereplicate_2:
	""" remove the duplicates inside the single file """
	input:
		config["output"] + "/06_concatenated_data/all_reads.fasta"
	output:
		config["output"] + "/07_derep_data/unique_reads.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing redundant reads for all reads"
	log:
		config["output"] + "/logs/07_dereplicate/all_reads.txt"
	shell:
		"vsearch --derep_fulllength {input} --output {output} --sizein --sizeout &>> {log}"

rule denoising:
	input:
		config["output"] + "/07_derep_data/unique_reads.fasta"
	output:
		config["output"] + "/08_ASVs/ASVs.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Generating ASVs"
	log:
		config["output"] + "/logs/08_ASVs/all_reads.txt"
	shell:
		"vsearch --cluster_unoise {input} --centroids {output} --relabel ASV --sizein --sizeout &>>{log}"

rule remove_chimeras:
	input:
		config["output"] + "/08_ASVs/ASVs.fasta"
	output:
		config["output"] + "/09_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing Chimeras"
	log:
		config["output"] + "/logs/09_ASVs_nonchimeras/all_reads.txt"
	shell:
		"vsearch -uchime3_denovo {input} --nonchimeras {output} --sizein --xsize &>> {log}"

rule generate_community_table:
	input:
		search = config["output"] + "/06_concatenated_data/all_reads.fasta",
		db = config["output"] + "/09_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	output:
		community_table = config["output"] + "/10_community_table/community_table.txt",
	params:
		options = " ".join(config["community_table_options"]),
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Generating community table"
	log:
		config["output"] + "/logs/10_community_table/all_reads.txt"
	shell:
		"vsearch --usearch_global {input.search} --db {input.db} {params.options} --otutabout {output.community_table} &>> {log}"

rule generate_report:
	input:
		community_table = config["output"] + "/10_community_table/community_table.txt",
		custom_mqc_config = "multiqc_config.yaml"
	output:
		config["output"] + "/11_report/multiqc_report.html"
	conda:
		os.path.join(workflow.basedir, "envs/mb_multiqc.yaml")
	params:
		output_dir = config["output"] + "/11_report/",
		log_dir = config["output"] + "/logs"
	message:
		"Generating MultiQC report"
	log:
		config["output"] + "/logs/11_MultiQC/multiqc.txt"
	shell:
		"multiqc {params.log_dir} -o {params.output_dir} --config {input.custom_mqc_config}"

rule taxonomy:
	input:
		ASVs = config["output"] + "/09_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	output:
		base = config["output"] + "/12_taxonomy/taxonomy.txt",
		plot = config["output"] + "/12_taxonomy/taxonomy.krona.txt"
	params:
		direct_db_lst = config["direct_dbs"],
		hierarchical_db = config["hierarchical_db"],
		threshold = config["classification_threshold"],
		keep_results = "False"
	threads:
		config["threads"]
	message:
		"Starting Multilevel Taxonomic Classification"
	conda:
		os.path.join(workflow.basedir, "envs/mb_taxonomy.yaml")
	log:
		config["output"] + "/logs/12_taxonomy/taxonomy.log"
	shell:
		"python3 scripts/multilvl_taxonomic_classification.py -d {params.direct_db_lst} -z {input.ASVs} -t {params.threshold} -o {output.base} -n {threads} -p {params.hierarchical_db} -k {params.keep_results} -l {log}"

rule krona:
	input:
		config["output"] + "/12_taxonomy/taxonomy.krona.txt"
	output:
		config["output"] + "/12_taxonomy/krona_plot.html"
	conda:
		os.path.join(workflow.basedir, "envs/mb_krona.yaml")
	message:
		"Creating Krona Plot"
	shell:
		"ktImportText -q {input} -o {output}"

rule merge_tables:
	input:
		community_table = config["output"] + "/10_community_table/community_table.txt",
		tax_table = config["output"] + "/12_taxonomy/taxonomy.txt"
	output:
		merged_table = config["output"] + "/13_merged/community_and_tax_merged.txt"
	message:
		"Merging community and taxonomy Tables"
	run:
		community_table = pd.read_csv(input.community_table, sep='\t', index_col = 0)
		tax_table = pd.read_csv(input.tax_table, sep = '\t', index_col = 1).iloc[:,1:]
		merged_table = pd.merge(community_table, tax_table, left_index=True, right_index=True)
		merged_table.to_csv(output.merged_table, sep='\t')
