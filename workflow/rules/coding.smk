rule rename_headers_for_dnoise:
    """Rename sequence headers for DnoisE

    Input file to dnoise can have only one ; char in Fasta headers
    """
    input:
        "results/06_derep/unique_reads.fasta",
    output:
        "results/06_derep/unique_reads_rename.fasta",
    threads: 1
    log:
        "logs/06_derep/rename_headers_for_dnoise.log",
    shell:
        "sed 's/;sample/_sample/' {input} | sed 's/;ee/_ee/' > {output} 2> {log};"


rule denoising_dnoise:
    """Denoise with DnoisE algorithm - for protein-coding sequences

    Alpha parameter should be chosen empirically based on entropy ratio changes.
    Pipeline is run with default alpha value; diagnostic plots also generated
    for range of alpha values. User should then update the alpha and minsize
    values in config file to rerun if necessary.
    """
    input:
        "results/06_derep/unique_reads_rename.fasta",
    output:
        denoised="results/07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.fasta",
    conda:
        "../envs/mb_dnoise.yaml"
    params:
        prefix=lambda wildcards: f"results/07_ASVs/ASVs_{wildcards.alpha}",
        frame=config["coding"]["frame"],
    threads: 4
    message:
        "Denoising with DnoisE"
    log:
        "logs/07_ASVs/denoising_dnoise.{alpha}.log",
    shell:
        """
        dnoise --fasta_input {input} -y --alpha {wildcards.alpha} \
        -x {params.frame} --cores {threads} --fasta_output {params.prefix} &> {log};
        """


rule calc_entropy_dnoise:
    """Calculate entropy ratio of 2nd and 3rd codon positions
    """
    input:
        "{prefix}.fasta",
    output:
        "{prefix}_entropy_values.csv",
    conda:
        "../envs/mb_dnoise.yaml"
    threads: 4
    message:
        "Calculating entropy per codon position with DnoisE"
    log:
        "{prefix}_entropy_values.log",
    params:
        frame=config["coding"]["frame"],
    shell:
        """
        dnoise --fasta_input {input} -g -x {params.frame} \
        --cores {threads} --csv_output {wildcards.prefix} &> {log};
        """


rule plot_entropy_ratio_vs_alpha:
    """Plot entropy ratio for range of values of denoising parameter alpha

    Generate diagnostic plot for user to choose final value of alpha, if
    default value is not suitable.
    """
    input:
        expand(
            "results/07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d_entropy_values.csv",
            alpha=config["dnoise_opts"]["alpha_range"],
        ),
    output:
        "results/07_ASVs/entropy_ratio_denoising_plot.png",
    params:
        param_vals=",".join([str(i) for i in config["dnoise_opts"]["alpha_range"]]),
        param_name="alpha",
        inputs=lambda wildcards, input: ",".join(input),
    conda:
        "../envs/mb_dnoise.yaml"
    log:
        "logs/07_ASVs/plot_entropy_ratio_vs_alpha.log",
    script:
        "../scripts/plot_entropy_ratio.py"


rule plot_entropy_ratio_vs_minsize:
    """Plot entropy ratio for range of values of minimum cluster size

    Generate diagnostic plot for user to choose final value of minsize, if
    default value is not suitable.
    """
    input:
        expand(
            "results/07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.minsize_{minsize}_entropy_values.csv",
            minsize=config["dnoise_opts"]["minsize_range"],
            alpha=config["denoising"]["alpha"],
        ),
    output:
        "results/07_ASVs/entropy_ratio_minsize_plot.png",
    params:
        param_vals=",".join([str(i) for i in config["dnoise_opts"]["minsize_range"]]),
        param_name="minsize",
        inputs=lambda wildcards, input: ",".join(input),
    conda:
        "../envs/mb_dnoise.yaml"
    log:
        "logs/07_ASVs/plot_entropy_ratio_vs_minsize.log",
    script:
        "../scripts/plot_entropy_ratio.py"


rule rename_denoised_ASVs:
    input:
        expand(
            "results/07_ASVs/ASVs_{alpha}_Adcorr_denoised_ratio_d.minsize_{minsize}.fasta",
            alpha=config["denoising"]["alpha"],
            minsize=config["denoising"]["minsize"],
        ),
    output:
        "results/07_ASVs/ASVs_dnoise.fasta",
    conda:
        "../envs/mb_vsearch.yaml"
    log:
        "logs/07_ASVs/rename_denoise_ASVs.log",
    threads: 1
    shell:
        """
        vsearch --sortbysize {input} --output {output} --relabel ASV \
        --sizein --sizeout --fasta_width 0 &> {log};
        """


rule screen_pseudogenes:
    """Remove sequences that have in frame stops and/or fail HMM screen"""
    input:
        "results/07_ASVs/ASVs_{method}.fasta",
    output:
        screened="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.fasta",
        stats="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_stats.tsv",
        hmmsearch="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hmmsearch.out",
        hist_hmm="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_hmm.png",
        hist_spf="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_spf.png",
        hist_mins="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_mins.png",
        mqc_hmm="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_hmm_mqc.json",
        mqc_spf="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_spf_mqc.json",
        mqc_mins="results/08_ASVs_screened/ASVs_{method}.no_pseudogenes.screen_hist_mins_mqc.json",
    wildcard_constraints:
        method=r"[a-z]+",
    conda:
        "../envs/mb_pseudogenes.yaml"
    log:
        "logs/07_ASVs/screen_pseudogenes.{method}.log",
    threads: 1
    params:
        hmm=config["coding"]["hmm"],
        code=config["coding"]["code"],
    shell:
        """
        pytransaln --input {input} --hmm {params.hmm} \
        --code {params.code} --out_hmmsearch {output.hmmsearch} \
        stats --out_screened {output.screened} --out_hist_hmm {output.hist_hmm} \
        --out_mqc_hmm {output.mqc_hmm} --out_mqc_spf {output.mqc_spf} \
        --out_mqc_mins {output.mqc_mins} \
        --out_stats {output.stats} --out_hist_spf {output.hist_spf} \
        --out_hist_mins {output.hist_mins} &> {log};
        """
