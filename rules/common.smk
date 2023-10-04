rule dereplicate_2:
	""" Remove the duplicates inside the single file """
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
		"logs/06_dereplicate/all_reads.txt"
	shell:
		"vsearch --threads {threads} --sizein --sizeout --fasta_width 0 --derep_fulllength {input} --output {output} &>> {log}"

rule denoising_unoise:
	"""Denoise with Unoise algorithm - for non-coding sequences

	Unoise algorithm as implemented in Vsearch. User can specify alpha and
	minsize parameters in the config file.
	"""
	input:
		"06_derep_data/unique_reads.fasta"
	output:
		"07_ASVs/ASVs.fasta"
	params:
		alpha=config["denoise_alpha"],
		minsize=config["denoise_minsize"],
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Generating ASVs with vsearch cluster_unoise"
	log:
		"logs/07_ASVs/all_reads.txt"
	shell:
		"""
		vsearch --cluster_unoise {input} \
		--threads {threads} --sizein --sizeout --fasta_width 0 \
		--minsize {params.minsize} --unoise_alpha {params.alpha} \
		--centroids {output} --relabel ASV &>>{log}
		"""

rule fasta_minsize:
	"""Filter sequence clusters by minimum size based on size= key in header

	Sequence header should contain '=' separated key-value pairs as expected
	from Vsearch and DnoisE input/outputs. The key 'size' is parsed to get
	cluster size, e.g. ">sequence_id;size=1234"
	"""
	input:
		"{prefix}.fasta"
	output:
		"{prefix}.minsize_{minsize}.fasta"
	conda: "../envs/mb_vsearch.yaml"
	threads: 1
	shell:
		"""
		vsearch --sizein --sizeout --fasta_width 0 --minsize {wildcards.minsize} --sortbysize {input} --output {output}
		"""

rule remove_chimeras:
	input:
		"07_ASVs/ASVs.fasta"
	output:
		"08_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Removing Chimeras"
	log:
		"logs/08_ASVs_nonchimeras/all_reads.txt"
	shell:
		"""
		vsearch --threads {threads} --sizein --sizeout --fasta_width 0 \
		-uchime3_denovo {input} --nonchimeras {output} {params.options} &>> {log}"
		"""

rule generate_community_table:
	input:
		search = "05_concatenated_data/all_reads.fasta",
		db = "08_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	output:
		community_table = "09_community_table/community_table.txt",
		community_table_biom = "09_community_table/community_table.biom",
	params:
		options = " ".join(config["community_table_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Generating community table"
	log:
		"logs/09_community_table/all_reads.txt"
	shell:
		"""
		vsearch --threads {threads} --sizein --sizeout {params.options} \
		--usearch_global {input.search} --db {input.db} \
		--otutabout {output.community_table} --biomout {output.community_table_biom} &>> {log}
		"""

rule taxonomy:
	input:
		ASVs = "08_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	output:
		base = "10_taxonomy/taxonomy.txt",
		plot = "10_taxonomy/taxonomy.krona.txt",
		stat_table_mqc = "10_taxonomy/stats_mqc.csv"
	params:
		keep_results = True,
		hierarchical_threshold = config["hierarchical_threshold"],
		script_path = os.path.join(workflow.basedir, "scripts/multilvl_taxonomic_classification.py")
	threads: workflow.cores
	message:
		"Starting Multilevel Taxonomic Classification"
	conda:
		"../envs/mb_taxonomy.yaml"
	log:
		"logs/10_taxonomy/taxonomy.log"
	script: "{params.script_path}"

rule krona:
	input:
		"10_taxonomy/taxonomy.krona.txt"
	output:
		"10_taxonomy/krona_plot.html"
	conda:
		"../envs/mb_krona.yaml"
	message:
		"Creating Krona Plot"
	shell:
		"ktImportText -q {input} -o {output}"

rule merge_tables:
	input:
		community_table = "09_community_table/community_table.txt",
		tax_table = "10_taxonomy/taxonomy.txt"
	output:
		merged_table = "11_merged/community_and_tax_merged.txt"
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
		community_table = "09_community_table/community_table.txt",
		custom_mqc_config = os.path.join(workflow.basedir, "multiqc_config.yaml"),
		stat_table_mqc = "10_taxonomy/stats_mqc.csv"
	output:
		"12_report/multiqc_report.html"
	conda:
		"../envs/mb_multiqc.yaml"
	params:
		output_dir = "12_report/",
		log_dir = "logs"
	message:
		"Generating MultiQC report"
	log:
		"logs/12_MultiQC/multiqc.txt"
	shell:
		"""
		multiqc {params.log_dir} {input.stat_table_mqc} -o {params.output_dir} \
		--config {input.custom_mqc_config} &> {log}
		"""

# vim: set noexpandtab:
