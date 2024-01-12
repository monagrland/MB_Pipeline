rule concat_libs_per_sample:
    """Concatenate multiple sequencing files from the same sample"""
    input:
        fw=lambda wildcards: reads_df[reads_df["sample"] == wildcards.sample]["fwd"],
        rv=lambda wildcards: reads_df[reads_df["sample"] == wildcards.sample]["rev"],
    output:
        fw=temp("results/01_trimmed/{sample}.concat.R1.fastq.gz"),  # TODO get file extension from input
        rv=temp("results/01_trimmed/{sample}.concat.R2.fastq.gz"),
    shell:
        """
        cat {input.fw} > {output.fw};
        cat {input.rv} > {output.rv};
        """


rule cutadapt:
    """Remove adapter sequences from the reads"""
    input:
        fw="results/01_trimmed/{sample}.concat.R1.fastq.gz",
        rv="results/01_trimmed/{sample}.concat.R2.fastq.gz",
    output:
        fw=temp("results/01_trimmed/{sample}.trim.R1.fastq.gz"),
        rv=temp("results/01_trimmed/{sample}.trim.R2.fastq.gz"),
    params:
        adapter_5p=config["adapter_trimming_options"]["5p"],
        adapter_3p=config["adapter_trimming_options"]["3p"],
        min_overlap=config["adapter_trimming_options"]["min_overlap"],
        others_common=config["adapter_trimming_options"]["others_common"],
    conda:
        "../envs/mb_cutadapt.yaml"
    threads: 1
    log:
        "logs/01_cutadapt/{sample}.txt",
    shell:
        """
        cutadapt --cores {threads} \
        -g {params.adapter_5p} -G {params.adapter_3p} -O {params.min_overlap} \
        {params.others_common} \
        -o {output.fw} -p {output.rv} {input.fw} {input.rv} &>>  {log}
        """


rule merge:
    """Merge paired end reads to a single file"""
    input:
        fw="results/01_trimmed/{sample}.trim.R1.fastq.gz",
        rv="results/01_trimmed/{sample}.trim.R2.fastq.gz",
    output:
        temp("results/02_merged/{sample}.merged.fastq.gz"),
    params:
        options=" ".join(config["merge_options"]),
    conda:
        "../envs/mb_vsearch.yaml"
    threads: 1
    log:
        "logs/02_merging/{sample}.txt",
    shell:
        """
        vsearch --fastq_mergepairs {input.fw} --reverse {input.rv} \
        --fastqout {output} {params.options} --relabel {wildcards.sample}_ \
        --label_suffix \;sample={wildcards.sample} --threads {threads} &>> {log}
        """


rule merge_mqc:
    """MultiQC summary of merge statistics"""
    input:
        expand("logs/02_merging/{sample}.txt", sample=samples),
    output:
        "logs/vsearch_fastq_mergepairs._mqc.json",
    params:
        script_path=os.path.join(workflow.basedir, "scripts/vsearch_logs_multiqc.py"),
    shell:
        """
        python {params.script_path} --format fastq_mergepairs --files {input} > {output}
        """


rule quality_filter:
    """Filter reads by quality scores"""
    input:
        "results/02_merged/{sample}.merged.fastq.gz",
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
