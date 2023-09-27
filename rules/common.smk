rule dereplicate_2:
	""" Remove the duplicates inside the single file """
	input:
		"05_concatenated_data/all_reads.fasta"
	output:
		"06_derep_data/unique_reads.fasta"
	params:
		options = " ".join(config["derep2_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Removing redundant reads for all reads"
	log:
		"logs/06_dereplicate/all_reads.txt"
	shell:
		"vsearch --threads {threads} --derep_fulllength {input} --output {output} {params.options} &>> {log}"

rule denoising_unoise:
	input:
		"06_derep_data/unique_reads.fasta"
	output:
		"07_ASVs/ASVs.fasta"
	params:
		options = " ".join(config["denoise_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Generating ASVs with vsearch cluster_unoise"
	log:
		"logs/07_ASVs/all_reads.txt"
	shell:
		"vsearch --threads {threads} --cluster_unoise {input} --centroids {output} --relabel ASV {params.options} &>>{log}"

rule rename_headers_for_dnoise:
	# input file to dnoise can have only one ; char in Fasta headers
	input:
		"06_derep_data/unique_reads.fasta"
	output:
		"06_derep_data/unique_reads_rename.fasta"
	threads: 1
	shell:
		"sed 's/;sample/_sample/' {input} | sed 's/;ee/_ee/' > {output};"

rule denoising_dnoise:
	input:
		"06_derep_data/unique_reads_rename.fasta"
	output:
		denoised="07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.fasta"
	conda:
		"../envs/mb_dnoise.yaml"
	params:
		prefix=lambda wildcards: f"07_ASVs/ASVs_{wildcards.alpha}",
		frame=3, # TODO move to config file
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
	input:
		"07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.fasta"
	output:
		"07_ASVs/ASVs_{alpha}_entropy_values.csv"
	conda:
		"../envs/mb_dnoise.yaml"
	params:
		prefix=lambda wildcards: f"07_ASVs/ASVs_{wildcards.alpha}",
	threads: 4
	message:
		"Calculating entropy per codon position with DnoisE"
	log:
		"logs/calc_entropy_dnoise.{alpha}.log"
	shell:
		"""
		dnoise --fasta_input {input} -g --cores {threads} --csv_output {params.prefix} &> {log};
		"""

# TODO Calculate entropy at different minimum abundances

rule remove_chimeras:
	input:
		"07_ASVs/ASVs.fasta"
	output:
		"08_ASVs_nonchimeras/ASVs_nonchimeras.fasta"
	params:
		options = " ".join(config["chimera_check_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads: workflow.cores
	message:
		"Removing Chimeras"
	log:
		"logs/08_ASVs_nonchimeras/all_reads.txt"
	shell:
		"vsearch --threads {threads} -uchime3_denovo {input} --nonchimeras {output} {params.options} &>> {log}"

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
		vsearch --threads {threads} --usearch_global {input.search} --db {input.db} \
		{params.options} --otutabout {output.community_table} --biomout {output.community_table_biom} &>> {log}
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
