import os
import json
import pandas as pd
import csv
import gzip
import shutil


fw_files = glob_wildcards(config["directory"] + "/{prefix}_R1_{suffix}.gz")


rule all:
	input:
		config["output"] + "/12_report/multiqc_report.html",
		config["output"] + "/11_merged/community_and_tax_merged.txt",
		config["output"] + "/10_taxonomy/krona_plot.html",

rule cutadapt:
	""" Rule to remove the Adapter Sequences from the reads """
	input:
		input_fw = config["directory"] + "/{prefix}_R1_{suffix}.gz",
		input_rv = config["directory"] + "/{prefix}_R2_{suffix}.gz"
	output:
		output_fw = config["output"] + "/01_trimmed_data/{prefix}_R1_{suffix}.gz",
		output_rv = config["output"] + "/01_trimmed_data/{prefix}_R2_{suffix}.gz"
	params:
		options = " ".join(config["adapter_trimming_options"]),
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

rule merge:
	""" Rule to merge paired end reads to a single file"""
	input:
		input_fw = config["output"] + "/01_trimmed_data/{prefix}_R1_{suffix}.gz",
		input_rv = config["output"] + "/01_trimmed_data/{prefix}_R2_{suffix}.gz"
	output:
		config["output"] + "/02_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	params:
		options = " ".join(config["merge_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
		basename = "{prefix}"
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Merging paired end reads for {params.filename_fw} and {params.filename_rv}"
	log:
		config["output"] + "/logs/02_merging/{prefix}_{suffix}.txt"
	shell:
		"vsearch --fastq_mergepairs {input.input_fw} --reverse {input.input_rv} --fastqout {output} {params.options} --relabel {params.basename}_ --label_suffix \;sample={params.basename} &>> {log}"

rule quality_filter:
	""" Rule for quality filtering """
	input:
		config["output"] + "/02_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	output:
		config["output"] + "/03_filtered_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	params:
		options = " ".join(config["filter_options"]),
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Executing Quality Filtering"
	log:
		config["output"] + "/logs/03_quality_filtering/{prefix}_{suffix}.txt"
	shell:
		"vsearch --fastq_filter {input} {params.options} --fastaout {output}"

rule dereplicate:
	""" rule to remove duplicates """
	input:
		config["output"] + "/03_filtered_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	output:
		config["output"] + "/04_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	params:
		filename = "{prefix}_{suffix}.fasta",
		options = " ".join(config["derep1_options"])
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing redundant reads for {params.filename}"
	log:
		config["output"] + "/logs/04_dereplicate/{prefix}_" + os.path.splitext("{suffix}")[0] + ".txt"
	shell:
		"vsearch --derep_fulllength {input} --output {output} {params.options} &>> {log}"

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
	params:
		options = " ".join(config["derep2_options"])
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing redundant reads for all reads"
	log:
		config["output"] + "/logs/06_dereplicate/all_reads.txt"
	shell:
		"vsearch --derep_fulllength {input} --output {output} {params.options} &>> {log}"

rule denoising:
	input:
		config["output"] + "/06_derep_data/unique_reads.fasta"
	output:
		config["output"] + "/07_ASVs/ASVs.fasta"
	params:
		options = " ".join(config["denoise_options"])
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Generating ASVs"
	log:
		config["output"] + "/logs/07_ASVs/all_reads.txt"
	shell:
		"vsearch --cluster_unoise {input} --centroids {output} --relabel ASV {params.options} &>>{log}"

rule remove_chimeras:
	input:
		config["output"] + "/07_ASVs/ASVs.fasta"
	output:
		config["output"] + "/08_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	params:
		options = " ".join(config["chimera_check_options"])
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Removing Chimeras"
	log:
		config["output"] + "/logs/08_ASVs_nonchimeras/all_reads.txt"
	shell:
		"vsearch -uchime3_denovo {input} --nonchimeras {output} {params.options} &>> {log}"

rule generate_community_table:
	input:
		search = config["output"] + "/05_concatenated_data/all_reads.fasta",
		db = config["output"] + "/08_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	output:
		community_table = config["output"] + "/09_community_table/community_table.txt",
	params:
		options = " ".join(config["community_table_options"]),
	conda:
		os.path.join(workflow.basedir, "envs/mb_vsearch.yaml")
	message:
		"Generating community table"
	log:
		config["output"] + "/logs/09_community_table/all_reads.txt"
	shell:
		"vsearch --usearch_global {input.search} --db {input.db} {params.options} --otutabout {output.community_table} &>> {log}"

rule taxonomy:
	input:
		ASVs = config["output"] + "/08_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	output:
		base = config["output"] + "/10_taxonomy/taxonomy.txt",
		plot = config["output"] + "/10_taxonomy/taxonomy.krona.txt"
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
		config["output"] + "/logs/10_taxonomy/taxonomy.log"
	shell:
		"python3 scripts/multilvl_taxonomic_classification.py -d {params.direct_db_lst} -z {input.ASVs} -t {params.threshold} -o {output.base} -n {threads} -p {params.hierarchical_db} -k {params.keep_results} -l {log}"

rule krona:
	input:
		config["output"] + "/10_taxonomy/taxonomy.krona.txt"
	output:
		config["output"] + "/10_taxonomy/krona_plot.html"
	conda:
		os.path.join(workflow.basedir, "envs/mb_krona.yaml")
	message:
		"Creating Krona Plot"
	shell:
		"ktImportText -q {input} -o {output}"

rule merge_tables:
	input:
		community_table = config["output"] + "/09_community_table/community_table.txt",
		tax_table = config["output"] + "/10_taxonomy/taxonomy.txt"
	output:
		merged_table = config["output"] + "/11_merged/community_and_tax_merged.txt"
	message:
		"Merging community and taxonomy Tables"
	run:
		community_table = pd.read_csv(input.community_table, sep='\t', index_col = 0)
		tax_table = pd.read_csv(input.tax_table, sep = '\t', index_col = 1).iloc[:,1:]
		merged_table = pd.merge(community_table, tax_table, left_index=True, right_index=True)
		merged_table.to_csv(output.merged_table, sep='\t')


rule generate_report:
	input:
		community_table = config["output"] + "/09_community_table/community_table.txt",
		custom_mqc_config = "multiqc_config.yaml"
	output:
		config["output"] + "/12_report/multiqc_report.html"
	conda:
		os.path.join(workflow.basedir, "envs/mb_multiqc.yaml")
	params:
		output_dir = config["output"] + "/12_report/",
		log_dir = config["output"] + "/logs"
	message:
		"Generating MultiQC report"
	log:
		config["output"] + "/logs/12_MultiQC/multiqc.txt"
	shell:
		"multiqc {params.log_dir} -o {params.output_dir} --config {input.custom_mqc_config}"
