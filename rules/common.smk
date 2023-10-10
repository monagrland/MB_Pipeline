wildcard_constraints:
		method=r"[a-z]+",
		screening=r"no_[a-z]+"

rule dereplicate_2:
	"""Remove exact duplicates in concatenated reads"""
	input:
		"05_concatenated_data/all_reads.fasta"
	output:
		"06_derep_data/unique_reads.fasta"
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Removing redundant reads for all reads"
	log:
		"logs/06_derep_data/dereplicate_2.log"
	shell:
		"""
		vsearch --derep_fulllength {input} --output {output} \
		--threads {threads} --sizein --sizeout --fasta_width 0 &>> {log}
		"""

rule denoising_unoise:
	"""Denoise with UNOISE algorithm - for non-coding sequences

	UNOISE algorithm as implemented in VSEARCH. User can specify alpha and
	minsize parameters in the config file.
	"""
	input:
		"06_derep_data/unique_reads.fasta"
	output:
		"07_ASVs/ASVs_unoise.fasta"
	params:
		alpha=config['denoising']['alpha'],
		minsize=config['denoising']["minsize"],
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Generating ASVs with vsearch cluster_unoise"
	log:
		"logs/07_ASVs/denoising_unoise.log"
	shell:
		"""
		vsearch --cluster_unoise {input} --centroids {output} --relabel ASV \
		--minsize {params.minsize} --unoise_alpha {params.alpha} \
		--threads {threads} --sizein --sizeout --fasta_width 0 &>>{log}
		"""

rule fasta_minsize:
	"""Filter sequence clusters by minimum size based on size= key in header

	Sequence header should contain '=' separated key-value pairs as expected
	from VSEARCH and DnoisE input/outputs. The key 'size' is parsed to get
	cluster size, e.g. ">sequence_id;size=1234"
	"""
	input:
		"{prefix}.fasta"
	output:
		"{prefix}.minsize_{minsize}.fasta"
	conda: "../envs/mb_vsearch.yaml"
	threads: 1
	log: "logs/{prefix}.fasta_minsize.{minsize}.log"
	shell:
		"""
		vsearch --minsize {wildcards.minsize} --sortbysize {input} --output {output} \
		--sizein --sizeout --fasta_width 0 &>> log;
		"""

rule remove_chimeras:
	"""Remove chimeras with UCHIME algorithm

	De-novo UNOISE algorithm as implemented in VSEARCH removes PCR chimera
	artefacts. If input sequences are protein-coding, the alternative
	pseudogene removal and HMM screening step should also remove PCR chimeras.
	"""
	input:
		"07_ASVs/ASVs_{method}.fasta"
	output:
		"08_ASVs_screened/ASVs_{method}.no_chimeras.fasta"
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Removing Chimeras"
	log:
		"logs/08_ASVs_screened/remove_chimeras.{method}.log"
	shell:
		"""
		vsearch -uchime3_denovo {input} --nonchimeras {output} \
		--threads {threads} --sizein --sizeout --fasta_width 0 &>> {log}
		"""

rule generate_community_table:
	"""Map reads to ASVs and generate counts per ASV per sample"""
	input:
		search = "05_concatenated_data/all_reads.fasta",
		db = "08_ASVs_screened/ASVs_{method}.{screening}.fasta"
	output:
		community_table = "09_community_table/community_table.{method}.{screening}.txt",
		community_table_biom = "09_community_table/community_table.{method}.{screening}.biom",
	params:
		options = " ".join(config["community_table_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Generating community table"
	log:
		"logs/09_community_table/generate_community_table.{method}.{screening}.log"
	shell:
		"""
		vsearch --usearch_global {input.search} --db {input.db} \
		--otutabout {output.community_table} --biomout {output.community_table_biom} \
		--threads {threads} --sizein --sizeout {params.options} &>> {log}
		"""

rule taxonomy:
	"""Assign taxonomy to ASVs with two-step taxonomic classification method"""
	input:
		ASVs = "08_ASVs_screened/ASVs_{method}.{screening}.fasta"
	output:
		base = "10_taxonomy/taxonomy.{method}.{screening}.txt",
		krona = "10_taxonomy/taxonomy.{method}.{screening}.krona.txt",
		stats_mqc = "10_taxonomy/stats_mqc.{method}.{screening}.csv"
	params:
		keep_results = True,
		hierarchical_threshold = config["hierarchical_threshold"],
		script_path = os.path.join(workflow.basedir, "scripts/multilvl_taxonomic_classification.py")
	threads: workflow.cores
	message:
		"Starting Multilevel Taxonomic Classification"
	conda:
		"../envs/mb_taxonomy.yaml"
	log: # TODO how to redirect script print output?
		"logs/10_taxonomy/taxonomy.{method}.{screening}.log"
	script: "{params.script_path}"

rule krona:
	"""Render Krona plot of the taxonomic summary"""
	input:
		"10_taxonomy/taxonomy.{method}.{screening}.krona.txt"
	output:
		"10_taxonomy/krona_plot.{method}.{screening}.html"
	conda:
		"../envs/mb_krona.yaml"
	message:
		"Creating Krona Plot"
	shell:
		"ktImportText -q {input} -o {output}"

rule merge_tables:
	"""Combine ASV counts per sample with taxonomy"""
	input:
		community_table = "09_community_table/community_table.{method}.{screening}.txt",
		tax_table = "10_taxonomy/taxonomy.{method}.{screening}.txt",
	output:
		merged_table = "11_merged/community_and_tax_merged.{method}.{screening}.txt"
	message:
		"Merging community and taxonomy Tables"
	run:
		community_table = pd.read_csv(input.community_table, sep='\t', index_col = 0)
		tax_table = pd.read_csv(input.tax_table, sep = '\t', index_col = 1).iloc[:,1:]
		ASV_with_size_lst = tax_table.index.tolist()
		ASV_sans_size_lst = []
		for ASV in ASV_with_size_lst:
		    ASV_sans_size_lst.append(ASV.split(";")[0])
		tax_table.index = ASV_sans_size_lst
		merged_table = pd.merge(community_table, tax_table, left_index=True, right_index=True)
		merged_table.to_csv(output.merged_table, sep='\t')

rule generate_report:
	input:
		community_table = "09_community_table/community_table.{method}.{screening}.txt",
		custom_mqc_config = os.path.join(workflow.basedir, "multiqc_config.yaml"),
		stats_mqc = "10_taxonomy/stats_mqc.{method}.{screening}.csv"
	output:
		"12_report/multiqc_report.{method}.{screening}.html"
	conda:
		"../envs/mb_multiqc.yaml"
	params:
		output_dir = lambda wildcards, output: os.path.dirname(output[0]),
		output_fn = lambda wildcards, output: os.path.basename(output[0]),
		log_dir = "logs"
	message:
		"Generating MultiQC report"
	log:
		"logs/12_report/generate_report.{method}.{screening}.log"
	shell:
		"""
		multiqc {input.stats_mqc} -n {params.output_fn} -o {params.output_dir} \
		--config {input.custom_mqc_config} {params.log_dir} --force &> {log}
		"""

# vim: set noexpandtab:
