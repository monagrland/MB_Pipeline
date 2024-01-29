wildcard_constraints:
    method=r"dnoise|unoise",
    screening=r"no_chimeras|no_pseudogenes",
    refdb=r"[^.]+",
    class_method=r"[^.]+",


rule quality_filter_mqc:
    input:
        expand("logs/03_quality_filtering/{sample}.txt", sample=samples),
    output:
        "logs/vsearch_fastq_filter._mqc.json",
    params:
        format="fastq_filter",
    log:
        "logs/03_quality_filtering/quality_filter_mqc.log",
    script:
        "../scripts/vsearch_logs_multiqc.py"


rule dereplicate:
    """Remove duplicates in each read file"""
    input:
        "results/03_filtered/{sample}.filtered.fasta",
    output:
        temp("results/04_derep/{sample}.derep.fasta"),
    conda:
        "../envs/mb_vsearch.yaml"
    threads: 1
    log:
        "logs/04_derep/{sample}.txt",
    shell:
        """
        vsearch --derep_fulllength {input} --output {output} \
        --threads {threads} --strand plus --sizeout --fasta_width 0 &>> {log}
        """


rule dereplicate_mqc:
    input:
        expand("logs/04_derep/{sample}.txt", sample=samples),
    output:
        "logs/vsearch_derep_fulllength._mqc.json",
    params:
        format="derep_fulllength",
    log:
        "logs/04_derep/dereplicate_mqc.log",
    script:
        "../scripts/vsearch_logs_multiqc.py"


rule concat_samples:
    """Concatenate reads from all samples into single file"""
    input:
        expand("results/04_derep/{sample}.derep.fasta", sample=samples),
    output:
        temp("results/05_concat/all_reads.fasta"),
    log:
        "logs/05_concat/all_reads.txt",
    shell:
        "cat {input} > {output} 2> {log}"


rule dereplicate_2:
    """Remove exact duplicates in concatenated reads"""
    input:
        "results/05_concat/all_reads.fasta",
    output:
        "results/06_derep/unique_reads.fasta",
    conda:
        "../envs/mb_vsearch.yaml"
    threads: workflow.cores
    message:
        "Removing redundant reads for all reads"
    log:
        "logs/06_derep/dereplicate_2.log",
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
        "results/06_derep/unique_reads.fasta",
    output:
        "results/07_ASVs/ASVs_unoise.fasta",
    params:
        alpha=config["denoising"]["alpha"],
        minsize=config["denoising"]["minsize"],
    conda:
        "../envs/mb_vsearch.yaml"
    threads: workflow.cores
    message:
        "Generating ASVs with vsearch cluster_unoise"
    log:
        "logs/07_ASVs/denoising_unoise.log",
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
        "{prefix}.fasta",
    output:
        "{prefix}.minsize_{minsize}.fasta",
    conda:
        "../envs/mb_vsearch.yaml"
    threads: 1
    log:
        "logs/{prefix}.fasta_minsize.{minsize}.log",
    shell:
        """
        vsearch --minsize {wildcards.minsize} --sortbysize {input} --output {output} \
        --sizein --sizeout --fasta_width 0 &>> {log};
        """


rule remove_chimeras:
    """Remove chimeras with UCHIME algorithm

    De-novo UNOISE algorithm as implemented in VSEARCH removes PCR chimera
    artefacts. If input sequences are protein-coding, the alternative
    pseudogene removal and HMM screening step should also remove PCR chimeras.
    """
    input:
        "results/07_ASVs/ASVs_{method}.fasta",
    output:
        "results/08_ASVs_screened/ASVs_{method}.no_chimeras.fasta",
    conda:
        "../envs/mb_vsearch.yaml"
    threads: workflow.cores
    message:
        "Removing Chimeras"
    log:
        "logs/08_ASVs_screened/remove_chimeras.{method}.log",
    shell:
        """
        vsearch -uchime3_denovo {input} --nonchimeras {output} \
        --threads {threads} --sizein --sizeout --fasta_width 0 &>> {log}
        """


rule generate_community_table:
    """Map reads to ASVs and generate counts per ASV per sample"""
    input:
        search="results/05_concat/all_reads.fasta",
        db="results/08_ASVs_screened/ASVs_{method}.{screening}.fasta",
    output:
        community_table="results/09_community_table/community_table.{method}.{screening}.txt",
        community_table_biom="results/09_community_table/community_table.{method}.{screening}.biom",
    params:
        options=" ".join(config["community_table_options"]),
    conda:
        "../envs/mb_vsearch.yaml"
    threads: workflow.cores
    message:
        "Generating community table"
    log:
        "logs/09_community_table/generate_community_table.{method}.{screening}.log",
    shell:  # --sizeout doesn't have effect here; size= doesn't appear in the otutab or biom files
        """
        vsearch --usearch_global {input.search} --db {input.db} \
        --otutabout {output.community_table} --biomout {output.community_table_biom} \
        --threads {threads} --sizein --sizeout {params.options} &>> {log}
        """


rule strip_size_info:
    """Remove abundance information from Fasta headers

    The ; and = characters used in key=value pairs cause problems with Newick
    tree files. Strip the ;size= abundance information from Fasta headers to
    keep only ASV identifiers.
    """
    input:
        "{prefix}.fasta",
    output:
        "{prefix}.nosize.fasta",
    conda:
        "../envs/mb_vsearch.yaml"
    log:
        "logs/{prefix}.strip_size_info.log",
    threads: 1
    shell:
        """
        vsearch --sortbysize {input} -xsize --output {output} &>> {log};
        """


rule align_screened_asvs:
    input:
        "results/08_ASVs_screened/ASVs_{method}.{screening}.nosize.fasta",
    output:
        "results/13_phylogeny/ASVs_{method}.{screening}.mafft",
    conda:
        "../envs/mb_phylogeny.yaml"
    log:
        "logs/13_phylogeny/align_screened_asvs.{method}.{screening}.log",
    threads: 8
    shell:
        """
        mafft --thread {threads} {input} > {output} 2>> {log};
        """


rule fasttree_screened_asvs:
    input:
        "results/13_phylogeny/ASVs_{method}.{screening}.mafft",
    output:
        "results/13_phylogeny/ASVs_{method}.{screening}.fasttree.treefile",
    conda:
        "../envs/mb_phylogeny.yaml"
    log:
        "logs/13_phylogeny/align_screened_asvs.{method}.{screening}.fasttree.log",
    threads: 1
    shell:
        """
        fasttree -nt -gtr -gamma {input} > {output} 2> {log};
        """


rule iqtree_screened_asvs:
    input:
        "results/13_phylogeny/ASVs_{method}.{screening}.mafft",
    output:
        "results/13_phylogeny/ASVs_{method}.{screening}.treefile",
    conda:
        "../envs/mb_phylogeny.yaml"
    log:
        "logs/13_phylogeny/align_screened_asvs.{method}.{screening}.iqtree.log",
    params:
        prefix="results/13_phylogeny/ASVs_{method}.{screening}.iqtree",
    threads: 8
    shell:
        """
        iqtree -s {input} -m TEST --fast --prefix {params.prefix} -redo -T AUTO --threads-max {threads} &>> {log}
        """


rule faiths_pd:
    input:
        tree="results/13_phylogeny/ASVs_{method}.{screening}.{treeprog}.treefile",
        biom="results/09_community_table/community_table.{method}.{screening}.biom",
    output:
        "results/13_phylogeny/ASVs_{method}.{screening}.{treeprog}.faiths_pd.tsv",
    conda:
        "../envs/mb_phylogeny.yaml"
    log:
        "logs/13_phylogeny/faiths_pd.{method}.{screening}.{treeprog}.log",
    script:
        "../scripts/faiths_pd.py"


rule taxonomy:
    """Assign taxonomy to ASVs with two-step taxonomic classification method"""
    input:
        ASVs="results/08_ASVs_screened/ASVs_{method}.{screening}.fasta",
        direct_dbs=lambda wildcards: config["dbpaths"]["multilvl"]["direct_dbs"][
            wildcards.refdb
        ],
        hierarchical_db=lambda wildcards: config["dbpaths"]["multilvl"][
            "hierarchical_db"
        ][wildcards.refdb],
    output:
        base="results/10_taxonomy/taxonomy.{method}.{screening}.multilvl.{refdb}.txt",
        krona="results/10_taxonomy/taxonomy.{method}.{screening}.multilvl.{refdb}.krona.txt",
        stats_mqc="results/10_taxonomy/stats_mqc.{method}.{screening}.multilvl.{refdb}.csv",
    params:
        keep_results=True,
        hierarchical_threshold=config["hierarchical_threshold"],
    threads: workflow.cores
    message:
        "Starting Multilevel Taxonomic Classification"
    conda:
        "../envs/mb_taxonomy.yaml"
    log:  # TODO how to redirect script print output?
        "logs/10_taxonomy/taxonomy.{method}.{screening}.multilvl.{refdb}.log",
    script:
        "../scripts/multilvl_taxonomic_classification.py"


rule sintax:
    """Assign taxonomy to ASVs with SINTAX implemented in Vsearch"""
    input:
        ASVs="results/08_ASVs_screened/ASVs_{method}.{screening}.fasta",
        db=lambda wildcards: config["dbpaths"]["sintax"][wildcards.refdb],
    output:
        sintax="results/10_taxonomy/taxonomy.{method}.{screening}.sintax.{refdb}.tsv",
    threads: workflow.cores
    conda:
        "../envs/mb_taxonomy.yaml"
    log:
        "logs/10_taxonomy/taxonomy.{method}.{screening}.sintax.{refdb}.log",
    shell:
        """
        vsearch --threads {threads} --db {input.db} --sintax {input.ASVs} --tabbedout {output} &> {log}
        """


rule sintax_to_krona:
    """Convert SINTAX output to Krona input table"""
    input:
        "results/10_taxonomy/taxonomy.{method}.{screening}.sintax.{refdb}.tsv",
    output:
        krona="results/10_taxonomy/taxonomy.{method}.{screening}.sintax.{refdb}.krona.txt",
        tax_table="results/10_taxonomy/taxonomy.{method}.{screening}.sintax.{refdb}.txt",
    conda:
        "../envs/mb_krona.yaml"
    log:
        "logs/10_taxonomy/sintax_to_krona.{method}.{screening}.sintax.{refdb}.log",
    params:
        cutoff=config["hierarchical_threshold"],
    script:
        "../scripts/sintax_output_krona.py"


rule krona_render:
    """Render Krona plot of the taxonomic summary"""
    input:
        "results/10_taxonomy/taxonomy.{method}.{screening}.{class_method}.{refdb}.krona.txt",
    output:
        "results/10_taxonomy/krona_plot.{method}.{screening}.{class_method}.{refdb}.html",
    conda:
        "../envs/mb_krona.yaml"
    message:
        "Creating Krona Plot"
    log:
        "logs/10_taxonomy/krona.{method}.{screening}.{class_method}.{refdb}.log",
    shell:
        "ktImportText -q {input} -o {output} &> {log}"


rule merge_tables:
    """Combine ASV counts per sample with taxonomy"""
    input:
        community_table="results/09_community_table/community_table.{method}.{screening}.txt",
        tax_table="results/10_taxonomy/taxonomy.{method}.{screening}.{class_method}.{refdb}.txt",
    output:
        "results/11_merged/community_and_tax_merged.{method}.{screening}.{class_method}.{refdb}.txt",
    log:
        "logs/11_merged/merge_tables.{method}.{screening}.{class_method}.{refdb}.log",
    script:
        "../scripts/merge_tables.py"


def mqc_files(paired, coding):
    out = [
        "logs/vsearch_fastq_filter._mqc.json",
        "logs/vsearch_derep_fulllength._mqc.json",
    ]
    if paired:
        out.append("logs/vsearch_fastq_mergepairs._mqc.json")
    if coding:
        out.extend(
            [
                "results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_hmm_mqc.json",
                "results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_spf_mqc.json",
                "results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_mins_mqc.json",
            ]
        )
    if "multilvl" in config["class_method"]:
        out.append(
            "results/10_taxonomy/stats_mqc.{method}.{screening}.multilvl.{refdb}.csv"
        )
    return out


rule generate_report:
    input:
        community_table="results/09_community_table/community_table.{method}.{screening}.txt",
        custom_mqc_config="config/multiqc_config.yaml",  # TODO
        mqc_files=mqc_files(config["paired"], config["protein_coding"]),
    output:
        "results/12_report/multiqc_report.{method}.{screening}.{refdb}.html",
    conda:
        "../envs/mb_multiqc.yaml"
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output[0]),
        output_fn=lambda wildcards, output: os.path.basename(output[0]),
        log_dir="logs",
    message:
        "Generating MultiQC report"
    log:
        "logs/12_report/generate_report.{method}.{screening}.{refdb}.log",
    shell:
        """
        multiqc {input.mqc_files} \
        -n {params.output_fn} -o {params.output_dir} \
        --config {input.custom_mqc_config} {params.log_dir} --force &> {log}
        """
