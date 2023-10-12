rule cutadapt:
	"""Remove the Adapter Sequences from the reads"""
	input:
		input_fw = os.path.join(config["input"], "{basename}.gz"),
	output:
		temp(output_fw = "01_trimmed_data/{basename}.gz"),
	params:
		options = " ".join(config["adapter_trimming_options"]),
	conda:
		"../envs/mb_cutadapt.yaml"
	threads: 1
	log:
		"logs/01_cutadapt/{basename}.txt"
	shell:
		"""
		cutadapt --cores {threads} {params.options} -o {output.output_fw} \
		{input.input_fw} &>>  {log}
		"""

rule relabel:
	"""Relabel fastq headers"""
	input:
		"01_trimmed_data/{basename}.gz",
	output:
		temp("02_relabeled/{basename}")
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
	"""Filter reads by quality scores"""
	input:
		"02_relabeled/{basename}"
	output:
		temp("03_filtered_data/{basename}.fasta")
	params:
		options = " ".join(config["filter_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Executing Quality Filtering"
	log:
		"logs/03_quality_filtering/{basename}.txt"
	shell:
		"""
		vsearch --fastq_filter {input} {params.options} --fastaout {output} \
		--threads {threads} &>> {log}
		"""

rule dereplicate:
	"""Remove duplicates and relabel the samples"""
	input:
		"03_filtered_data/{basename}.fasta"
	output:
		temp("04_derep_data/{basename}.fasta")
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	log:
		"logs/04_dereplicate/{basename}.txt"
	shell:
		"""
		vsearch --derep_fulllength {input} --output {output} \
		--threads {threads} --strand plus --sizeout --fasta_width 0 &>> {log}
		"""

rule concatenate:
	"""Concatenate all reads into single file"""
	input:
		expand(
			"04_derep_data/{basename}.fasta",
			zip,
			basename = files_single.basename,
		)
	output:
		temp("05_concatenated_data/all_reads.fasta")
	params:
		derep_dir = "04_derep_data/"
	message:
		"Concatenating all reads"
	log:
		"logs/05_concatenate/all_reads.txt"
	shell:
		"cat {params.derep_dir}*.fasta > {output} 2> {log}"

# vim: set noexpandtab:
