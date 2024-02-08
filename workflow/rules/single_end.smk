rule concat_libs_per_sample:
    """Concatenate multiple sequencing files from the same sample"""
    input:
        fw=lambda wildcards: reads_df[reads_df["sample"] == wildcards.sample]["fwd"],
    output:
        temp("results/01_trimmed/{sample}.concat.fastq.gz"),  # TODO get file extension from input
    shell:
        """
        cat {input.fw} > {output};
        """


rule cutadapt:
    """Remove adapter sequences from the reads"""
    input:
        "results/01_trimmed/{sample}.concat.fastq.gz",  # TODO get file extension from input
    output:
        temp("results/01_trimmed/{sample}.trim.fastq.gz"),
    params:
        adapter_5p=config["adapter_trimming_options"]["5p"],
        min_overlap=config["adapter_trimming_options"]["min_overlap"],
        others_common=config["adapter_trimming_options"]["others_common"],
    conda:
        "../envs/mb_cutadapt.yaml"
    threads: 1
    log:
        "logs/01_cutadapt/{sample}.txt",
    shell:
        """
        cutadapt --cores {threads} -g {params.adapter_5p} -O {params.min_overlap} \
        {params.others_common} \
        -o {output} {input} &>> {log}
        """


rule relabel:
    """Relabel Fastq headers"""
    input:
        "results/01_trimmed/{sample}.trim.fastq.gz",
    output:
        temp("results/02_relabeled/{sample}.fastq.gz"),
    params:
        script_path=os.path.join(workflow.basedir, "scripts/relabel.py"),
    threads: 1
    shell:
        """
        python3 {params.script_path} -i {input} -o {output}
        """


rule quality_filter_single:
    """Filter reads by quality scores"""
    input:
        "results/02_relabeled/{sample}.fastq.gz",
    output:
        temp("results/03_filtered/{sample}.filtered.fasta"),
    params:
        options=" ".join(config["filter_options"]),
    conda:
        "../envs/mb_vsearch.yaml"
    threads: 1
    log:
        "logs/03_quality_filtering/{sample}.txt",
    shell:
        """
        vsearch --fastq_filter {input} {params.options} --fastaout {output} \
        --threads {threads} &>> {log}
        """
