rule cutadapt:
	""" Rule to remove the Adapter Sequences from the reads """
	input:
		input_fw = os.path.join(config["directory"], "{prefix}_R1_{suffix}.gz"),
		input_rv = os.path.join(config["directory"], "{prefix}_R2_{suffix}.gz")
	output:
		output_fw = temp(os.path.join(config["output"], "01_trimmed_data/{prefix}_R1_{suffix}.gz")),
		output_rv = temp(os.path.join(config["output"], "01_trimmed_data/{prefix}_R2_{suffix}.gz"))
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
		os.path.join(config["output"], "logs/01_cutadapt/{prefix}_{suffix}.txt")
	shell:
		"""
		cutadapt --cores {threads} {params.options} -o {output.output_fw} -p {output.output_rv} \
		{input.input_fw} {input.input_rv} &>>  {log}
		"""

rule merge:
	""" Rule to merge paired end reads to a single file"""
	input:
		input_fw = os.path.join(config["output"], "01_trimmed_data/{prefix}_R1_{suffix}.gz"),
		input_rv = os.path.join(config["output"], "01_trimmed_data/{prefix}_R2_{suffix}.gz")
	output:
		temp(os.path.join(config["output"], "02_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fastq"))
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
		os.path.join(config["output"], "logs/02_merging/{prefix}_{suffix}.txt")
	shell:
		"""
		vsearch --threads {threads} \
		--fastq_mergepairs {input.input_fw} --reverse {input.input_rv} \
		--fastqout {output} {params.options} --relabel {params.basename}_ \
		--label_suffix \;sample={params.basename} &>> {log}
		"""

rule quality_filter:
	""" Rule for quality filtering """
	input:
		os.path.join(config["output"], "02_merged_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fastq")
	output:
		temp(os.path.join(config["output"], "03_filtered_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"))
	params:
		options = " ".join(config["filter_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Executing Quality Filtering"
	log:
		os.path.join(config["output"], "logs/03_quality_filtering/{prefix}_{suffix}.txt")
	shell:
		"vsearch --threads {threads} --fastq_filter {input} {params.options} --fastaout {output} &>> {log}"

rule dereplicate:
	""" Rule to remove duplicates """
	input:
		os.path.join(config["output"], "03_filtered_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta")
	output:
		temp(os.path.join(config["output"], "04_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"))
	params:
		filename = "{prefix}_{suffix}.fasta",
		options = " ".join(config["derep1_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads: 1
	message:
		"Removing redundant reads for {params.filename}"
	log:
		os.path.join(config["output"], "logs/04_dereplicate/{prefix}_" + os.path.splitext("{suffix}")[0] + ".txt")
	shell:
		"vsearch --threads {threads} --derep_fulllength {input} --output {output} {params.options} &>> {log}"

rule concatenate:
	""" Rule to concatenate all files into one """
	input:
		expand(
			os.path.join(
				config["output"],
				"04_derep_data/{prefix}_" + os.path.splitext("{suffix}")[0] + "_merged.fasta"
			),
			zip,
			prefix = fw_files.prefix,
			suffix = fw_files.suffix
		)
	output:
		temp(os.path.join(config["output"], "05_concatenated_data/all_reads.fasta"))
	params:
		derep_dir = os.path.join(config["output"], "04_derep_data/")
	message:
		"Concatenating all reads"
	log:
		os.path.join(config["output"], "logs/05_concatenate/all_reads.txt")
	shell:
		"cat {params.derep_dir}*.fasta > {output} 2> {log}"
