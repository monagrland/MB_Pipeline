rule cutadapt:
	"""Remove the Adapter Sequences from the reads"""
	input:
		input_fw = os.path.join(config["input"], "{prefix}_R1_{suffix}.gz"),
		input_rv = os.path.join(config["input"], "{prefix}_R2_{suffix}.gz")
	output:
		output_fw = temp("01_trimmed_data/{prefix}_R1_{suffix}.gz"),
		output_rv = temp("01_trimmed_data/{prefix}_R2_{suffix}.gz")
	params:
		options = " ".join(config["adapter_trimming_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
	conda:
		"../envs/mb_cutadapt.yaml"
	threads: 1
	message:
		"Executing adapter trimming for {params.filename_fw} and {params.filename_rv}"
	log:
		"logs/01_cutadapt/{prefix}_{suffix}.txt"
	shell:
		"""
		cutadapt --cores {threads} {params.options} -o {output.output_fw} \
		-p {output.output_rv} {input.input_fw} {input.input_rv} &>>  {log}
		"""

rule merge:
	"""Merge paired end reads to a single file"""
	input:
		input_fw = "01_trimmed_data/{prefix}_R1_{suffix}.gz",
		input_rv = "01_trimmed_data/{prefix}_R2_{suffix}.gz"
	output:
		temp("02_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fastq")
	params:
		options = " ".join(config["merge_options"]),
		filename_fw = "{prefix}_R1_{suffix}.gz",
		filename_rv = "{prefix}_R2_{suffix}.gz",
		basename = "{prefix}"
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Merging paired end reads for {params.filename_fw} and {params.filename_rv}"
	log:
		"logs/02_merging/{prefix}_{suffix}.txt"
	shell:
		"""
		vsearch --fastq_mergepairs {input.input_fw} --reverse {input.input_rv} \
		--fastqout {output} {params.options} --relabel {params.basename}_ \
		--label_suffix \;sample={params.basename} --threads {threads} &>> {log}
		"""

rule quality_filter:
	"""Filter reads by quality scores"""
	input:
		"02_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fastq"
	output:
		temp("03_filtered_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta")
	params:
		options = " ".join(config["filter_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Executing Quality Filtering"
	log:
		"logs/03_quality_filtering/{prefix}_{suffix}.txt"
	shell:
		"vsearch --fastq_filter {input} {params.options} --fastaout {output} --threads {threads} &>> {log}"

rule dereplicate:
	"""Remove duplicates in each read file"""
	input:
		"03_filtered_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
	output:
		temp("04_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta")
	params:
		filename = "{prefix}_{suffix}.fasta",
		options = " ".join(config["derep1_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Removing redundant reads for {params.filename}"
	log:
		"logs/04_dereplicate/{prefix}_" + os.path.splitext("{suffix}")[0] + ".txt"
	shell:
		"vsearch --threads {threads} --derep_fulllength {input} --output {output} {params.options} &>> {log}"

rule concatenate:
	"""Concatenate all reads into single file"""
	input:
		expand(
			"04_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta",
			zip,
			prefix = fw_files.prefix,
			suffix = fw_files.suffix
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
