rule rename_headers_for_dnoise:
	"""Rename sequence headers for DnoisE

	Input file to dnoise can have only one ; char in Fasta headers
	"""
	input:
		"06_derep_data/unique_reads.fasta"
	output:
		"06_derep_data/unique_reads_rename.fasta"
	threads: 1
	shell:
		"sed 's/;sample/_sample/' {input} | sed 's/;ee/_ee/' > {output};"

rule denoising_dnoise:
	"""Denoise with DnoisE algorithm - for protein-coding sequences

	Alpha parameter should be chosen empirically based on entropy ratio changes.
	Pipeline is run with default alpha value; diagnostic plots also generated
	for range of alpha values. User should then update the alpha and minsize
	values in config file to rerun if necessary.
	"""
	input:
		"06_derep_data/unique_reads_rename.fasta"
	output:
		denoised="07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.fasta"
	conda:
		"../envs/mb_dnoise.yaml"
	params:
		prefix=lambda wildcards: f"07_ASVs/ASVs_{wildcards.alpha}",
		frame=config["reading_frame_start"],
	threads: 4
	message:
		"Denoising with DnoisE"
	log:
		"logs/denoising_dnoise.{alpha}.log"
	shell:
		"""
		dnoise --fasta_input {input} -y --alpha {wildcards.alpha} -x {params.frame} --cores {threads} --fasta_output {params.prefix} &> {log};
		"""

rule calc_entropy_dnoise:
	"""Calculate entropy ratio of 2nd and 3rd codon positions
	"""
	input:
		"{prefix}.fasta"
	output:
		"{prefix}_entropy_values.csv"
	conda:
		"../envs/mb_dnoise.yaml"
	threads: 4
	message:
		"Calculating entropy per codon position with DnoisE"
	log:
		"{prefix}_entropy_values.log"
	params:
		frame=config["reading_frame_start"]
	shell:
		"""
		dnoise --fasta_input {input} -g -x {params.frame} --cores {threads} --csv_output {wildcards.prefix} &> {log};
		"""

rule plot_entropy_ratio_vs_alpha:
	"""Plot entropy ratio for range of values of denoising parameter alpha

	Generate diagnostic plot for user to choose final value of alpha, if
	default value is not suitable.
	"""
	input:
		expand("07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d_entropy_values.csv", alpha=config["alpha_range"])
	output:
		"diagnostics/entropy_ratio_denoising_plot.png"
	params:
		alphas=",".join([str(i) for i in config["alpha_range"]]),
		inputs=lambda wildcards, input: ",".join(input),
		script_path = os.path.join(workflow.basedir, "scripts/plot_entropy_ratio.py"),
	conda:
		"../envs/mb_dnoise.yaml"
	script: "{params.script_path}"

rule plot_entropy_ratio_vs_minsize:
	"""Plot entropy ratio for range of values of minimum cluster size

	Generate diagnostic plot for user to choose final value of minsize, if
	default value is not suitable.
	"""
	input:
		expand(
			"07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.minsize_{minsize}_entropy_values.csv",
			minsize=config['minsize_range'],
			alpha=config["denoise_alpha"]
		)
	output:
		"diagnostics/entropy_ratio_minsize_plot.png"
	params: # TODO update script arguments
		alphas=",".join([str(i) for i in config['minsize_range']]),
		inputs=lambda wildcards, input: ",".join(input),
		script_path = os.path.join(workflow.basedir, "scripts/plot_entropy_ratio.py"),
	conda:
		"../envs/mb_dnoise.yaml"
	script: "{params.script_path}"

rule rename_denoised_ASVs:
	input:
		expand("07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.minsize_{minsize}.fasta", alpha=config["denoise_alpha"], minsize=config["denoise_minsize"])
	output:
		"07_ASVs/ASVs_dnoise.fasta"
	conda: "../envs/mb_vsearch.yaml"
	threads: 1
	shell:
		"""
		vsearch --sizein --sizeout --fasta_width 0 --sortbysize {input} --output {output} --relabel ASV;
		"""

# vim: set noexpandtab:
