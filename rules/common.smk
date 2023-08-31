rule dereplicate_2:
	""" Remove the duplicates inside the single file """
	input:
		os.path.join(config["output"], "05_concatenated_data/all_reads.fasta")
	output:
		os.path.join(config["output"], "06_derep_data/unique_reads.fasta")
	params:
		options = " ".join(config["derep2_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads:
		config["threads"]
	message:
		"Removing redundant reads for all reads"
	log:
		os.path.join(config["output"], "logs/06_dereplicate/all_reads.txt")
	shell:
		"vsearch --threads {threads} --derep_fulllength {input} --output {output} {params.options} &>> {log}"

rule denoising:
	input:
		os.path.join(config["output"], "06_derep_data/unique_reads.fasta")
	output:
		os.path.join(config["output"], "07_ASVs/ASVs.fasta")
	params:
		options = " ".join(config["denoise_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads:
		config["threads"]
	message:
		"Generating ASVs"
	log:
		os.path.join(config["output"], "logs/07_ASVs/all_reads.txt")
	shell:
		"vsearch --threads {threads} --cluster_unoise {input} --centroids {output} --relabel ASV {params.options} &>>{log}"

rule remove_chimeras:
	input:
		os.path.join(config["output"], "07_ASVs/ASVs.fasta")
	output:
		os.path.join(config["output"], "08_ASVs_nonchimeras/ASVs_nonchimeras.fasta")
	params:
		options = " ".join(config["chimera_check_options"])
	conda:
		"../envs/mb_vsearch.yaml"
	threads:
		config["threads"]
	message:
		"Removing Chimeras"
	log:
		os.path.join(config["output"], "logs/08_ASVs_nonchimeras/all_reads.txt")
	shell:
		"vsearch --threads {threads} -uchime3_denovo {input} --nonchimeras {output} {params.options} &>> {log}"

rule generate_community_table:
	input:
		search = os.path.join(config["output"], "05_concatenated_data/all_reads.fasta"),
		db = os.path.join(config["output"], "08_ASVs_nonchimeras/ASVs_nonchimeras.fasta")
	output:
		community_table = os.path.join(config["output"], "09_community_table/community_table.txt"),
		community_table_biom = os.path.join(config["output"], "09_community_table/community_table.biom"),
	params:
		options = " ".join(config["community_table_options"]),
	conda:
		"../envs/mb_vsearch.yaml"
	threads:
		config["threads"]
	message:
		"Generating community table"
	log:
		os.path.join(config["output"], "logs/09_community_table/all_reads.txt")
	shell:
		"""
		vsearch --threads {threads} --usearch_global {input.search} --db {input.db} \
		{params.options} --otutabout {output.community_table} --biomout {output.community_table_biom} &>> {log}
		"""

rule taxonomy:
	input:
		ASVs = os.path.join(config["output"], "08_ASVs_nonchimeras/ASVs_nonchimeras.fasta")
	output:
		base = os.path.join(config["output"], "10_taxonomy/taxonomy.txt"),
		plot = os.path.join(config["output"], "10_taxonomy/taxonomy.krona.txt"),
		stat_table_mqc = os.path.join(config["output"], "10_taxonomy/stats_mqc.csv")
	params:
		keep_results = False,
		hierarchical_threshold = config["hierarchical_threshold"],
		script_path = os.path.join(workflow.basedir, "scripts/multilvl_taxonomic_classification.py")
	threads:
		config["threads"]
	message:
		"Starting Multilevel Taxonomic Classification"
	conda:
		"../envs/mb_taxonomy.yaml"
	log:
		os.path.join(config["output"], "logs/10_taxonomy/taxonomy.log")
	script: "{params.script_path}"

rule krona:
	input:
		os.path.join(config["output"], "10_taxonomy/taxonomy.krona.txt")
	output:
		os.path.join(config["output"], "10_taxonomy/krona_plot.html")
	conda:
		"../envs/mb_krona.yaml"
	message:
		"Creating Krona Plot"
	shell:
		"ktImportText -q {input} -o {output}"

rule merge_tables:
	input:
		community_table = os.path.join(config["output"], "09_community_table/community_table.txt"),
		tax_table = os.path.join(config["output"], "10_taxonomy/taxonomy.txt")
	output:
		merged_table = os.path.join(config["output"], "11_merged/community_and_tax_merged.txt")
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
		community_table = os.path.join(config["output"], "09_community_table/community_table.txt"),
		custom_mqc_config = os.path.join(workflow.basedir, "multiqc_config.yaml"),
		stat_table_mqc = os.path.join(config["output"], "10_taxonomy/stats_mqc.csv")
	output:
		os.path.join(config["output"], "12_report/multiqc_report.html")
	conda:
		"../envs/mb_multiqc.yaml"
	params:
		output_dir = os.path.join(config["output"], "12_report/"),
		log_dir = os.path.join(config["output"], "logs")
	message:
		"Generating MultiQC report"
	log:
		os.path.join(config["output"], "logs/12_MultiQC/multiqc.txt")
	shell:
		"""
		multiqc {params.log_dir} {input.stat_table_mqc} -o {params.output_dir} \
		--config {input.custom_mqc_config} &> {log}
		"""
