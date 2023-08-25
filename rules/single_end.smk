rule cutadapt:
	""" Rule to remove the Adapter Sequences from the reads """
	input:
		input_fw = os.path.join(config["directory"], "{basename}.gz"),

	output:
		output_fw = os.path.join(config["output"], "01_trimmed_data/{basename}.gz"),

	params:
		options = " ".join(config["adapter_trimming_options"]),
		filename_fw = "{basename}.gz",
	conda:
		"../envs/mb_cutadapt.yaml"
	threads: 1
	message:
		"Executing adapter trimming for {params.filename_fw}"
	log:
		os.path.join(config["output"], "logs/01_cutadapt/{basename}.txt")
	shell:
		"""
		cutadapt --cores {threads} {params.options} -o {output.output_fw} \
		{input.input_fw} &>>  {log}
		"""

rule relabel:
	""" Rule for relabeling of the fastq headers """
	input:
		os.path.join(config["output"], "01_trimmed_data/{basename}.gz"),
	output:
		os.path.join(config["output"], "02_relabeled/{basename}")
	params:
		script_path = os.path.join(workflow.basedir, "scripts/relabel.py")
	threads: 1
	message:
		"Relabeling FASTQ Headers"
	shell:
		"""
		python3 {params.script_path} -i {input} -o {output}
		"""

rule quality_filter_single:
	""" Rule for quality filtering """
	input:
		os.path.join(config["output"], "02_relabeled/{basename}")
	output:
		os.path.join(config["output"], "03_filtered_data/{basename}.fasta")
	params:
		options = " ".join(config["filter_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Executing Quality Filtering"
	log:
		os.path.join(config["output"], "logs/03_quality_filtering/{basename}.txt")
	shell:
		"""
		vsearch --threads {threads} --fastq_filter {input} {params.options} --fastaout {output} &>> {log}
		"""

rule dereplicate:
	""" Rule to remove duplicates and relabel the samples"""
	input:
		os.path.join(config["output"], "03_filtered_data/{basename}.fasta")
	output:
		os.path.join(config["output"], "04_derep_data/{basename}.fasta")
	params:
		filename = "{basename}.fasta",
		options = " ".join(config["derep1_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Removing redundant reads for {params.filename}"
	log:
		os.path.join(config["output"], "logs/04_dereplicate/{basename}.txt")
	shell:
		"vsearch --threads {threads} --derep_fulllength {input} --output {output} {params.options} &>> {log}"

rule concatenate:
	""" Rule to concatenate all files into one """
	input:
		expand(
			os.path.join(config["output"], "04_derep_data/{basename}.fasta"),
			zip,
			basename = files_single.basename,
		)
	output:
		os.path.join(config["output"], "05_concatenated_data/all_reads.fasta")
	params:
		derep_dir = os.path.join(config["output"], "04_derep_data/")
	message:
		"Concatenating all reads"
	log:
		os.path.join(config["output"], "logs/05_concatenate/all_reads.txt")
	shell:
		"cat {params.derep_dir}*.fasta > {output} 2> {log}"
