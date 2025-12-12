# RNAModBench - A comprehensive nanopore RNA modification detection pipeline
# Version: 1.0.0

import os
import glob
from pathlib import Path

# Pipeline configuration
configfile: "config/config.yaml"

# Sample and tool definitions
SAMPLES = config["samples"]
TOOLS = config["tools"]
REFERENCE_DIR = config["reference_dir"]
DATA_DIR = config["data_dir"]
RESULTS_DIR = config["results_dir"]

# Wildcard constraints
wildcard_constraints:
    sample="[^/]+",
    tool="[^/]+"

# Target rule - run complete analysis
def get_all_targets(wildcards):
    targets = []
    for sample in SAMPLES:
        for tool in TOOLS:
            if tool in ["CHEUI", "DENA", "DRUMMER", "ELIGOS2", "m6Anet", "Nanocompore", 
                       "MINES", "Epinano_SVM", "Epinano_DiffErr", "Nanom6A", "xPore", "yanocomp", "NanoSPA"]:
                if tool == "Epinano_DiffErr":
                    # Epinano DiffErr needs control sample
                    for control in SAMPLES:
                        if control != sample:
                            targets.append(f"{RESULTS_DIR}/{tool}/{sample}_vs_{control}/{sample}_vs_{control}_Epinano_DiffErr_processed.txt")
                elif tool in ["DRUMMER", "xPore", "yanocomp"]:
                    # These tools need paired analysis
                    for control in SAMPLES:
                        if control != sample:
                            targets.append(f"{RESULTS_DIR}/{tool}/{sample}_vs_{control}/{sample}_vs_{control}_{tool}_processed.txt")
                else:
                    targets.append(f"{RESULTS_DIR}/{tool}/{sample}/{sample}_{tool}_processed.txt")
    
    # Add visualization targets
    targets.extend([
        f"{RESULTS_DIR}/summary/depth_plots.png",
        f"{RESULTS_DIR}/summary/guitar_plots_mrna.png",
        f"{RESULTS_DIR}/summary/guitar_plots_ncrna.png",
        f"{RESULTS_DIR}/summary/liftover_summary.tsv"
    ])
    
    return targets

rule all:
    input:
        get_all_targets

# Rule to create output directories
rule create_directories:
    output:
        directory(f"{RESULTS_DIR}/logs"),
        directory(f"{RESULTS_DIR}/temp")
    shell:
        "mkdir -p {RESULTS_DIR}/logs {RESULTS_DIR}/temp"

# Data preprocessing rules
rule basecalling:
    input:
        fast5_dir=f"{DATA_DIR}/{{sample}}/fast5"
    output:
        fastq=f"{RESULTS_DIR}/basecalling/{{sample}}_basecalled.fastq",
        summary=f"{RESULTS_DIR}/basecalling/{{sample}}_sequencing_summary.txt"
    params:
        config=config["guppy"]["config"],
        callers=config["guppy"]["num_callers"],
        threads=config["guppy"]["threads_per_caller"]
    threads: config["guppy"]["total_threads"]
    shell:
        """
        guppy_basecaller --input_path {input.fast5_dir} \
            --recursive \
            --save_path {RESULTS_DIR}/basecalling/{wildcards.sample} \
            --fast5_out \
            --config {params.config} \
            --num_callers {params.callers} \
            --cpu_threads_per_caller {params.threads} \
            --device auto
        
        cat {RESULTS_DIR}/basecalling/{wildcards.sample}/pass/*.fastq > {output.fastq}
        cp {RESULTS_DIR}/basecalling/{wildcards.sample}/sequencing_summary.txt {output.summary}
        """

rule alignment_transcriptome:
    input:
        fastq=f"{RESULTS_DIR}/basecalling/{{sample}}_basecalled.fastq",
        transcriptome=f"{REFERENCE_DIR}/transcriptome.fa"
    output:
        sam=f"{RESULTS_DIR}/alignment/{{sample}}_transcriptome.sam",
        bam=f"{RESULTS_DIR}/alignment/{{sample}}_transcriptome.bam",
        bai=f"{RESULTS_DIR}/alignment/{{sample}}_transcriptome.bam.bai"
    threads: config["alignment"]["threads"]
    shell:
        """
        minimap2 -ax map-ont --MD -t {threads} {input.transcriptome} {input.fastq} > {output.sam}
        
        samtools view -@ {threads} -bh -F 2324 {output.sam} | \
            samtools sort -@ {threads} -o {output.bam}
        
        samtools index -@ {threads} {output.bam}
        """

rule alignment_genome:
    input:
        fastq=f"{RESULTS_DIR}/basecalling/{{sample}}_basecalled.fastq",
        genome=f"{REFERENCE_DIR}/genome.fa"
    output:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}_genome.bam",
        sorted_bam=f"{RESULTS_DIR}/alignment/{{sample}}_genome_sorted.bam",
        bai=f"{RESULTS_DIR}/alignment/{{sample}}_genome_sorted.bam.bai"
    threads: config["alignment"]["threads"]
    shell:
        """
        minimap2 -ax splice -k14 -t {threads} {input.genome} {input.fastq} | \
            samtools view -hSb | samtools sort -@ {threads} -o {output.bam}
        
        samtools view {output.bam} -bh -t {input.genome}.fai -F 2308 | \
            samtools sort -@ {threads} -o {output.sorted_bam}
        
        samtools index -@ {threads} {output.sorted_bam}
        """

rule nanopolish_index:
    input:
        fast5_dir=f"{DATA_DIR}/{{sample}}/fast5",
        fastq=f"{RESULTS_DIR}/basecalling/{{sample}}_basecalled.fastq"
    output:
        index_done=touch(f"{RESULTS_DIR}/nanopolish/{{sample}}_indexed.txt")
    shell:
        """
        nanopolish index -d {input.fast5_dir} {input.fastq}
        """

rule nanopolish_eventalign:
    input:
        fastq=f"{RESULTS_DIR}/basecalling/{{sample}}_basecalled.fastq",
        bam=f"{RESULTS_DIR}/alignment/{{sample}}_transcriptome.bam",
        genome=f"{REFERENCE_DIR}/transcriptome.fa",
        summary=f"{RESULTS_DIR}/basecalling/{{sample}}_sequencing_summary.txt",
        index_done=f"{RESULTS_DIR}/nanopolish/{{sample}}_indexed.txt"
    output:
        eventalign=f"{RESULTS_DIR}/nanopolish/{{sample}}_eventalign.txt"
    threads: config["nanopolish"]["threads"]
    shell:
        """
        nanopolish eventalign --reads {input.fastq} \
            --bam {input.bam} \
            --genome {input.genome} \
            --scale-events \
            --summary {input.summary} \
            --signal-index \
            --threads {threads} > {output.eventalign}
        """

# CHEUI analysis
rule cheui_preprocess:
    input:
        eventalign=f"{RESULTS_DIR}/nanopolish/{{sample}}_eventalign.txt"
    output:
        preprocessed=f"{RESULTS_DIR}/CHEUI/{{sample}}/preprocess_m6A/{sample}_eventalign_signals+IDS.p"
    params:
        kmer_model=config["cheui"]["kmer_model"],
        threads=config["cheui"]["threads"]
    conda:
        "envs/cheui.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/CHEUI/{wildcards.sample}/preprocess_m6A
        python -m sklearnex {workflow.basedir}/scripts/CHEUI_preprocess_m6A.py \
            -i {input.eventalign} \
            -m {params.kmer_model} \
            -n {params.threads} \
            -o {RESULTS_DIR}/CHEUI/{wildcards.sample}/preprocess_m6A
        """

rule cheui_predict_model1:
    input:
        preprocessed=f"{RESULTS_DIR}/CHEUI/{{sample}}/preprocess_m6A/{sample}_eventalign_signals+IDS.p"
    output:
        predictions=f"{RESULTS_DIR}/CHEUI/{{sample}}/read_level_m6A_predictions1.txt"
    params:
        model=config["cheui"]["model1"],
        label=f"{wildcards.sample}_m6A"
    conda:
        "envs/cheui.yaml"
    shell:
        """
        CUDA_VISIBLE_DEVICES=0 python -m sklearnex {workflow.basedir}/scripts/CHEUI_predict_model1.py \
            -i {input.preprocessed} \
            -m {params.model} \
            -o {output.predictions} \
            -l {params.label}
        """

rule cheui_sort_predictions:
    input:
        predictions=f"{RESULTS_DIR}/CHEUI/{{sample}}/read_level_m6A_predictions1.txt"
    output:
        sorted=f"{RESULTS_DIR}/CHEUI/{{sample}}/read_level_m6A_sorted.txt"
    shell:
        """
        sort -k1 --parallel=40 {input.predictions} > {output.sorted}
        """

rule cheui_predict_model2:
    input:
        sorted=f"{RESULTS_DIR}/CHEUI/{{sample}}/read_level_m6A_sorted.txt"
    output:
        final=f"{RESULTS_DIR}/CHEUI/{{sample}}/{sample}_CHEUI_raw.txt"
    params:
        model=config["cheui"]["model2"]
    conda:
        "envs/cheui.yaml"
    shell:
        """
        CUDA_VISIBLE_DEVICES=0 python -m sklearnex {workflow.basedir}/scripts/CHEUI_predict_model2.py \
            -i {input.sorted} \
            -m {params.model} \
            -o {output.final}
        """

rule cheui_postprocess:
    input:
        raw=f"{RESULTS_DIR}/CHEUI/{{sample}}/{sample}_CHEUI_raw.txt"
    output:
        processed=f"{RESULTS_DIR}/CHEUI/{{sample}}/{sample}_CHEUI_processed.txt"
    params:
        prob_threshold=config["cheui"]["prob_threshold"],
        ratio_threshold=config["cheui"]["ratio_threshold"]
    script:
        "scripts/postprocess_cheui.py"

# ELIGOS2 analysis
rule eligos2_solo:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}_genome_sorted.bam",
        bed=f"{REFERENCE_DIR}/genes.bed",
        ref=f"{REFERENCE_DIR}/genome.fa"
    output:
        result=f"{RESULTS_DIR}/ELIGOS2_solo/{{sample}}/{sample}_ELIGOS2_solo_raw.txt"
    params:
        prefix=wildcards.sample,
        outdir=f"{RESULTS_DIR}/ELIGOS2_solo/{wildcards.sample}",
        threads=config["eligos2"]["threads"],
        max_depth=config["eligos2"]["max_depth"],
        min_depth=config["eligos2"]["min_depth"]
    conda:
        "envs/eligos2.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        eligos2 rna_mod \
            -i {input.bam} \
            -reg {input.bed} \
            -ref {input.ref} \
            -p {params.prefix} \
            -o {params.outdir} \
            --max_depth {params.max_depth} \
            --min_depth {params.min_depth} \
            -t {params.threads}
        """

rule eligos2_postprocess:
    input:
        raw=f"{RESULTS_DIR}/ELIGOS2_solo/{{sample}}/{sample}_ELIGOS2_solo_raw.txt"
    output:
        processed=f"{RESULTS_DIR}/ELIGOS2_solo/{{sample}}/{sample}_ELIGOS2_solo_processed.txt"
    params:
        padj_threshold=config["eligos2"]["padj_threshold"],
        oddr_threshold=config["eligos2"]["oddr_threshold"]
    script:
        "scripts/postprocess_eligos2.py"

# m6Anet analysis
rule m6anet_dataprep:
    input:
        eventalign=f"{RESULTS_DIR}/nanopolish/{{sample}}_eventalign.txt"
    output:
        dataprep_dir=directory(f"{RESULTS_DIR}/m6Anet/{{sample}}_dataprep")
    params:
        processes=config["m6anet"]["dataprep_threads"],
        readcount_max=config["m6anet"]["readcount_max"]
    conda:
        "envs/m6anet.yaml"
    shell:
        """
        m6anet dataprep \
            --eventalign {input.eventalign} \
            --out_dir {output.dataprep_dir} \
            --n_processes {params.processes} \
            --readcount_max {params.readcount_max}
        """

rule m6anet_inference:
    input:
        dataprep_dir=f"{RESULTS_DIR}/m6Anet/{{sample}}_dataprep"
    output:
        inference_dir=directory(f"{RESULTS_DIR}/m6Anet/{{sample}}_inference")
    params:
        processes=config["m6anet"]["inference_threads"]
    conda:
        "envs/m6anet.yaml"
    shell:
        """
        m6anet inference \
            --input_dir {input.dataprep_dir} \
            --out_dir {output.inference_dir} \
            --n_processes {params.processes}
        """

rule m6anet_postprocess:
    input:
        inference_dir=f"{RESULTS_DIR}/m6Anet/{{sample}}_inference"
    output:
        processed=f"{RESULTS_DIR}/m6Anet/{{sample}}/{sample}_m6Anet_processed.txt"
    params:
        prob_threshold=config["m6anet"]["prob_threshold"],
        ratio_threshold=config["m6anet"]["ratio_threshold"]
    script:
        "scripts/postprocess_m6anet.py"

# Nanocompore analysis
rule nanocompore_eventalign_collapse:
    input:
        eventalign=f"{RESULTS_DIR}/nanopolish/{{sample}}_eventalign.txt"
    output:
        collapsed_dir=directory(f"{RESULTS_DIR}/Nanocompore/{{sample}}_collapse")
    params:
        threads=config["nanocompore"]["threads"]
    conda:
        "envs/nanocompore.yaml"
    shell:
        """
        nanocompore eventalign_collapse \
            -i {input.eventalign} \
            -o {output.collapsed_dir} \
            -t {params.threads}
        """

rule nanocompore_sampcomp:
    input:
        collapsed1=f"{RESULTS_DIR}/Nanocompore/{{sample}}_collapse/out_eventalign_collapse.tsv",
        collapsed2=f"{RESULTS_DIR}/Nanocompore/{{control}}_collapse/out_eventalign_collapse.tsv",
        fasta=f"{REFERENCE_DIR}/transcriptome.fa",
        bed=f"{REFERENCE_DIR}/genes.bed"
    output:
        result_dir=directory(f"{RESULTS_DIR}/Nanocompore/{{sample}}_vs_{{control}}")
    params:
        label1=wildcards.sample,
        label2=wildcards.control,
        min_coverage=config["nanocompore"]["min_coverage"],
        min_ref_length=config["nanocompore"]["min_ref_length"],
        threads=config["nanocompore"]["threads"]
    conda:
        "envs/nanocompore.yaml"
    shell:
        """
        nanocompore sampcomp \
            --file_list1 {input.collapsed1} \
            --file_list2 {input.collapsed2} \
            --label1 {params.label1} \
            --label2 {params.label2} \
            --fasta {input.fasta} \
            --bed {input.bed} \
            --outpath {output.result_dir} \
            --min_coverage {params.min_coverage} \
            --min_ref_length {params.min_ref_length} \
            --nthreads {params.threads} \
            --overwrite \
            --allow_warnings
        """

rule nanocompore_postprocess:
    input:
        result_dir=f"{RESULTS_DIR}/Nanocompore/{{sample}}_vs_{{control}}"
    output:
        processed=f"{RESULTS_DIR}/Nanocompore/{{sample}}_vs_{{control}}/{sample}_vs_{control}_Nanocompore_processed.txt"
    params:
        pvalue_threshold=config["nanocompore"]["pvalue_threshold"],
        lor_threshold=config["nanocompore"]["lor_threshold"]
    script:
        "scripts/postprocess_nanocompore.py"

# DENA analysis
rule dena_extract_pos:
    input:
        fasta=f"{REFERENCE_DIR}/transcriptome.fa"
    output:
        pos_file=f"{RESULTS_DIR}/DENA/{wildcards.sample}/{sample}_motif_positions.txt"
    params:
        motif=config["dena"]["motif"]
    conda:
        "envs/dena.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/DENA/{wildcards.sample}
        python -m sklearnex {workflow.basedir}/scripts/DENA_extract.py get_pos \
            --fasta {input.fasta} \
            --motif {params.motif} \
            --output {output.pos_file}
        """

rule dena_predict:
    input:
        fast5_dir=f"{DATA_DIR}/{{sample}}/fast5",
        bam=f"{RESULTS_DIR}/alignment/{{sample}}_transcriptome.bam",
        pos_file=f"{RESULTS_DIR}/DENA/{{sample}}/{sample}_motif_positions.txt"
    output:
        tmp_dir=directory(f"{RESULTS_DIR}/DENA/{{sample}}/tmp")
    params:
        corr_grp=config["dena"]["corr_grp"],
        windows=config["dena"]["windows"],
        processes=config["dena"]["processes"]
    conda:
        "envs/dena.yaml"
    shell:
        """
        python -m sklearnex {workflow.basedir}/scripts/DENA_extract.py predict \
            --fast5 {input.fast5_dir} \
            --corr_grp {params.corr_grp} \
            --bam {input.bam} \
            --sites {input.pos_file} \
            --label {wildcards.sample} \
            --windows {params.windows} \
            --processes {params.processes}
        """

rule dena_lstm_predict:
    input:
        tmp_dir=f"{RESULTS_DIR}/DENA/{{sample}}/tmp"
    output:
        result=f"{RESULTS_DIR}/DENA/{{sample}}/{sample}_DENA_raw.txt"
    params:
        model=config["dena"]["model"],
        prefix=wildcards.sample
    conda:
        "envs/dena.yaml"
    shell:
        """
        python -m sklearnex {workflow.basedir}/scripts/DENA_LSTM_predict.py \
            -i {input.tmp_dir} \
            -m {params.model} \
            -o {output.result} \
            -p {params.prefix} \
            -d
        """

rule dena_postprocess:
    input:
        raw=f"{RESULTS_DIR}/DENA/{{sample}}/{sample}_DENA_raw.txt"
    output:
        processed=f"{RESULTS_DIR}/DENA/{{sample}}/{sample}_DENA_processed.txt"
    params:
        ratio_threshold=config["dena"]["ratio_threshold"],
        coverage_threshold=config["dena"]["coverage_threshold"]
    script:
        "scripts/postprocess_dena.py"

# Epinano analysis
rule epinano_variants:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}_genome_sorted.bam",
        ref=f"{REFERENCE_DIR}/genome.fa"
    output:
        fwd_csv=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_fwd.per.site.csv",
        rev_csv=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_rev.per.site.csv"
    params:
        threads=config["epinano"]["threads"]
    conda:
        "envs/epinano.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/Epinano/{wildcards.sample}
        python {workflow.basedir}/scripts/Epinano_Variants.py \
            -c {params.threads} \
            -r {input.ref} \
            -b {input.bam}
        
        awk -F',' 'NR==1 || $5 >= 20' {wildcards.sample}_fwd.per.site.csv > {output.fwd_csv}
        awk -F',' 'NR==1 || $5 >= 20' {wildcards.sample}_rev.per.site.csv > {output.rev_csv}
        """

rule epinano_slide:
    input:
        fwd_csv=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_fwd.per.site.csv",
        rev_csv=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_rev.per.site.csv"
    output:
        fwd_5mer=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_fwd.5mer.csv",
        rev_5mer=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_rev.5mer.csv"
    conda:
        "envs/epinano.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/Slide_Variants.py {input.fwd_csv} 5
        python {workflow.basedir}/scripts/Slide_Variants.py {input.rev_csv} 5
        
        mv {wildcards.sample}_fwd.per.site.5mer.csv {output.fwd_5mer}
        mv {wildcards.sample}_rev.per.site.5mer.csv {output.rev_5mer}
        """

rule epinano_predict:
    input:
        fwd_5mer=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_fwd.5mer.csv",
        rev_5mer=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_rev.5mer.csv"
    output:
        fwd_pred=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_fwd.prediction.csv",
        rev_pred=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_rev.prediction.csv"
    params:
        model=config["epinano"]["model"],
        columns=config["epinano"]["columns"]
    conda:
        "envs/epinano.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/Epinano_Predict.py \
            --model {params.model} \
            --predict {input.fwd_5mer} \
            --columns {params.columns} \
            --out_prefix fwd_{wildcards.sample}
        
        python {workflow.basedir}/scripts/Epinano_Predict.py \
            --model {params.model} \
            --predict {input.rev_5mer} \
            --columns {params.columns} \
            --out_prefix rev_{wildcards.sample}
        """

rule epinano_merge:
    input:
        fwd_pred=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_fwd.prediction.csv",
        rev_pred=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_rev.prediction.csv"
    output:
        merged=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_Epinano_raw.txt"
    shell:
        """
        cat {input.fwd_pred} {input.rev_pred} > {output.merged}
        """

rule epinano_postprocess:
    input:
        raw=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_Epinano_raw.txt"
    output:
        processed=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_Epinano_processed.txt"
    params:
        delta_threshold=config["epinano"]["delta_threshold"]
    script:
        "scripts/postprocess_epinano.py"

# MINES analysis
rule tombo_resquiggle:
    input:
        fast5_dir=f"{DATA_DIR}/{{sample}}/fast5",
        ref=f"{REFERENCE_DIR}/transcriptome.fa"
    output:
        done=touch(f"{RESULTS_DIR}/Tombo/{{sample}}_resquiggled.txt")
    params:
        corr_grp=config["tombo"]["corr_grp"],
        basecall_grp=config["tombo"]["basecall_grp"],
        processes=config["tombo"]["processes"]
    conda:
        "envs/tombo.yaml"
    shell:
        """
        tombo resquiggle {input.fast5_dir} {input.ref} \
            --rna \
            --corrected-group {params.corr_grp} \
            --basecall-group {params.basecall_grp} \
            --overwrite \
            --processes {params.processes} \
            --fit-global-scale \
            --ignore-read-locks
        """

rule tombo_fraction:
    input:
        fast5_dir=f"{DATA_DIR}/{{sample}}/fast5",
        ref=f"{REFERENCE_DIR}/transcriptome.fa"
    output:
        fraction_wig=f"{RESULTS_DIR}/Tombo/{{sample}}/{sample}_fraction_modified_reads.plus.wig",
        coverage_bed=f"{RESULTS_DIR}/Tombo/{{sample}}/{sample}_coverage.plus.bedgraph"
    params:
        corr_grp=config["tombo"]["corr_grp"]
    conda:
        "envs/tombo.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/Tombo/{wildcards.sample}
        tombo detect_modifications \
            --fast5-basedirs {input.fast5_dir} \
            --statistics-file-basename {RESULTS_DIR}/Tombo/{wildcards.sample}/{wildcards.sample} \
            --rna \
            --corrected-group {params.corr_grp}
        """

rule mines_convert:
    input:
        fraction_wig=f"{RESULTS_DIR}/Tombo/{{sample}}/{sample}_fraction_modified_reads.plus.wig",
        coverage_bed=f"{RESULTS_DIR}/Tombo/{{sample}}/{sample}_coverage.plus.bedgraph"
    output:
        fraction_bed=f"{RESULTS_DIR}/MINES/{{sample}}/{sample}_fraction_modified_reads.plus.bed"
    shell:
        """
        mkdir -p {RESULTS_DIR}/MINES/{wildcards.sample}
        wig2bed < {input.fraction_wig} > {output.fraction_bed}
        """

rule mines_predict:
    input:
        fraction_bed=f"{RESULTS_DIR}/MINES/{{sample}}/{sample}_fraction_modified_reads.plus.bed",
        coverage_bed=f"{RESULTS_DIR}/Tombo/{{sample}}/{sample}_coverage.plus.bedgraph"
    output:
        mines_bed=f"{RESULTS_DIR}/MINES/{{sample}}/{sample}_MINES_raw.bed"
    params:
        ref=f"{REFERENCE_DIR}/transcriptome.fa",
        kmer_models=config["mines"]["kmer_models"]
    conda:
        "envs/mines.yaml"
    shell:
        """
        python {workflow.basedir}/scripts/MINES_cDNA.py \
            --fraction_modified {input.fraction_bed} \
            --coverage {input.coverage_bed} \
            --output {output.mines_bed} \
            --ref {params.ref} \
            --kmer_models {params.kmer_models}
        """

rule mines_postprocess:
    input:
        raw=f"{RESULTS_DIR}/MINES/{{sample}}/{sample}_MINES_raw.bed"
    output:
        processed=f"{RESULTS_DIR}/MINES/{{sample}}/{sample}_MINES_processed.txt"
    params:
        coverage_threshold=config["mines"]["coverage_threshold"],
        ratio_threshold=config["mines"]["ratio_threshold"]
    script:
        "scripts/postprocess_mines.py"

# xPore analysis
rule xpore_dataprep:
    input:
        eventalign=f"{RESULTS_DIR}/nanopolish/{{sample}}_eventalign.txt"
    output:
        dataprep_dir=directory(f"{RESULTS_DIR}/xPore/{{sample}}_dataprep")
    params:
        processes=config["xpore"]["dataprep_threads"],
        readcount_max=config["xpore"]["readcount_max"]
    conda:
        "envs/xpore.yaml"
    shell:
        """
        xpore dataprep \
            --eventalign {input.eventalign} \
            --out_dir {output.dataprep_dir} \
            --n_processes {params.processes} \
            --readcount_max {params.readcount_max}
        """

rule xpore_diffmod:
    input:
        dataprep1=f"{RESULTS_DIR}/xPore/{{sample}}_dataprep",
        dataprep2=f"{RESULTS_DIR}/xPore/{{control}}_dataprep",
        config=f"{RESULTS_DIR}/xPore/{{sample}}_vs_{{control}}_config.yml"
    output:
        diffmod_dir=directory(f"{RESULTS_DIR}/xPore/{{sample}}_vs_{{control}}_diffmod")
    conda:
        "envs/xpore.yaml"
    shell:
        """
        xpore diffmod --config {input.config} --n_processes 12
        """

rule xpore_postprocess:
    input:
        diffmod_dir=f"{RESULTS_DIR}/xPore/{{sample}}_vs_{{control}}_diffmod"
    output:
        processed=f"{RESULTS_DIR}/xPore/{{sample}}_vs_{{control}}/{sample}_vs_{control}_xPore_processed.txt"
    params:
        prob_threshold=config["xpore"]["prob_threshold"]
    script:
        "scripts/postprocess_xpore.py"

# yanocomp analysis
rule yanocomp_prep:
    input:
        eventalign=f"{RESULTS_DIR}/nanopolish/{{sample}}_eventalign.txt",
        gtf=f"{REFERENCE_DIR}/genes.gtf"
    output:
        hdf5=f"{RESULTS_DIR}/yanocomp/{{sample}}/{sample}.hdf5"
    params:
        processes=config["yanocomp"]["prep_threads"]
    conda:
        "envs/yanocomp.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/yanocomp/{wildcards.sample}
        yanocomp prep \
            -p {params.processes} \
            -e {input.eventalign} \
            -h {output.hdf5} \
            -g {input.gtf}
        """

rule yanocomp_gmmtest:
    input:
        control_hdf5=f"{RESULTS_DIR}/yanocomp/{{control}}/{control}.hdf5",
        treat_hdf5=f"{RESULTS_DIR}/yanocomp/{{sample}}/{sample}.hdf5"
    output:
        bed=f"{RESULTS_DIR}/yanocomp/{{sample}}_vs_{{control}}/{sample}_vs_{control}_yanocomp_raw.bed"
    params:
        processes=config["yanocomp"]["test_threads"],
        n_components=config["yanocomp"]["n_components"],
        fdr=config["yanocomp"]["fdr"]
    conda:
        "envs/yanocomp.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/yanocomp/{wildcards.sample}_vs_{wildcards.control}
        yanocomp gmmtest \
            -c {input.control_hdf5} \
            -t {input.treat_hdf5} \
            -p {params.processes} \
            -n {params.n_components} \
            -o {output.bed} \
            -f {params.fdr}
        """

rule yanocomp_postprocess:
    input:
        raw=f"{RESULTS_DIR}/yanocomp/{{sample}}_vs_{{control}}/{sample}_vs_{control}_yanocomp_raw.bed"
    output:
        processed=f"{RESULTS_DIR}/yanocomp/{{sample}}_vs_{{control}}/{sample}_vs_{control}_yanocomp_processed.txt"
    params:
        pvalue_threshold=config["yanocomp"]["pvalue_threshold"]
    script:
        "scripts/postprocess_yanocomp.py"

# NanoSPA analysis
rule nanospa:
    input:
        fast5_dir=f"{DATA_DIR}/{{sample}}/fast5",
        ref=f"{REFERENCE_DIR}/genome.fa"
    output:
        result_dir=directory(f"{RESULTS_DIR}/NanoSPA/{{sample}}")
    conda:
        "envs/nanospa.yaml"
    shell:
        """
        mkdir -p {output.result_dir}
        cd {output.result_dir}
        
        nanospa alignment -i {input.fast5_dir} -r {input.ref}
        nanospa remove_intron
        nanospa extract_features
        nanospa preprocess_m6A
        nanospa prediction_m6A
        nanospa prediction_psU
        """

rule nanospa_postprocess:
    input:
        result_dir=f"{RESULTS_DIR}/NanoSPA/{{sample}}"
    output:
        processed=f"{RESULTS_DIR}/NanoSPA/{{sample}}/{sample}_NanoSPA_processed.txt"
    params:
        prob_threshold=config["nanospa"]["prob_threshold"]
    script:
        "scripts/postprocess_nanospa.py"

# Utility rules
rule extract_5mer:
    input:
        bed=f"{{prefix}}.bed",
        fasta=f"{REFERENCE_DIR}/genome.fa"
    output:
        bed5mer=f"{{prefix}}_5mer.bed"
    script:
        "scripts/extract_5mer.py"

rule liftover:
    input:
        bed=f"{{prefix}}.bed",
        gtf=f"{REFERENCE_DIR}/genes.gtf"
    output:
        lifted=f"{{prefix}}_liftover.bed"
    params:
        r2d_tool=config["utilities"]["r2d_tool"]
    shell:
        """
        {params.r2d_tool} liftover -H -g {input.gtf} -i {input.bed} > {output.lifted}
        """

# Summary and comparison rules
rule generate_summary:
    input:
        expand("{RESULTS_DIR}/{tool}/{sample}/{sample}_{tool}_processed.txt", 
               RESULTS_DIR=RESULTS_DIR, 
               sample=SAMPLES, 
               tool=TOOLS)
    output:
        summary=f"{RESULTS_DIR}/summary/modification_summary.tsv",
        comparison=f"{RESULTS_DIR}/summary/tool_comparison.tsv"
    script:
        "scripts/generate_summary.py"

rule create_report:
    input:
        summary=f"{RESULTS_DIR}/summary/modification_summary.tsv",
        comparison=f"{RESULTS_DIR}/summary/tool_comparison.tsv"
    output:
        report=f"{RESULTS_DIR}/report/RNAModBench_report.html"
    script:
        "scripts/create_report.py"

# Depth analysis rules
rule samtools_depth:
    input:
        bam=f"{RESULTS_DIR}/alignment/{{sample}}_genome_sorted.bam"
    output:
        depth=f"{RESULTS_DIR}/depth/{{sample}}_depth.txt"
    shell:
        """
        mkdir -p {RESULTS_DIR}/depth
        samtools depth {input.bam} > {output.depth}
        """

rule generate_depth_plot:
    input:
        depth=f"{RESULTS_DIR}/depth/{{sample}}_depth.txt"
    output:
        plot=f"{RESULTS_DIR}/plots/{{sample}}_depth.png"
    params:
        sample=wildcards.sample
    conda:
        "envs/r.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/plots
        Rscript {workflow.basedir}/scripts/generate_depth_plots.R {input.depth} {output.plot} {params.sample}
        """

# DRUMMER analysis rules
rule drummer_postprocess:
    input:
        summary=f"{RESULTS_DIR}/DRUMMER/{{prefix}}/summary.txt"
    output:
        processed=f"{RESULTS_DIR}/DRUMMER/{{prefix}}/{wildcards.prefix}_DRUMMER_processed.txt"
    params:
        pvalue_threshold=config["drummer"]["pvalue_threshold"],
        frac_diff_threshold=config["drummer"]["frac_diff_threshold"]
    script:
        "scripts/postprocess_drummer.py"

# Epinano DiffErr analysis rules
rule epinano_differr:
    input:
        fwd1=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_fwd.per.site.csv",
        rev1=f"{RESULTS_DIR}/Epinano/{{sample}}/{sample}_rev.per.site.csv",
        fwd2=f"{RESULTS_DIR}/Epinano/{{control}}/{control}_fwd.per.site.csv",
        rev2=f"{RESULTS_DIR}/Epinano/{{control}}/{control}_rev.per.site.csv"
    output:
        merged=f"{RESULTS_DIR}/Epinano_DiffErr/{{sample}}_vs_{{control}}/{sample}_vs_{control}_Epinano_DiffErr_raw.txt"
    conda:
        "envs/epinano.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/Epinano_DiffErr/{wildcards.sample}_vs_{wildcards.control}
        cd {RESULTS_DIR}/Epinano_DiffErr/{wildcards.sample}_vs_{wildcards.control}
        
        # Run Epinano DiffErr for forward strand
        Rscript {workflow.basedir}/scripts/Epinano_DiffErr.R \
            -k {input.fwd1} \
            -w {input.fwd2} \
            -o {wildcards.sample}_vs_{wildcards.control}_fwd \
            -f sum_err
        
        # Run Epinano DiffErr for reverse strand
        Rscript {workflow.basedir}/scripts/Epinano_DiffErr.R \
            -k {input.rev1} \
            -w {input.rev2} \
            -o {wildcards.sample}_vs_{wildcards.control}_rev \
            -f sum_err
        
        # Merge results
        cat {wildcards.sample}_vs_{wildcards.control}_fwd.delta-sum_err.prediction.csv \
            {wildcards.sample}_vs_{wildcards.control}_rev.delta-sum_err.prediction.csv \
            > {output.merged}
        """

rule epinano_differr_postprocess:
    input:
        raw=f"{RESULTS_DIR}/Epinano_DiffErr/{{sample}}_vs_{{control}}/{sample}_vs_{control}_Epinano_DiffErr_raw.txt"
    output:
        processed=f"{RESULTS_DIR}/Epinano_DiffErr/{{sample}}_vs_{{control}}/{sample}_vs_{control}_Epinano_DiffErr_processed.txt"
    params:
        delta_threshold=config["epinano"]["delta_threshold"]
    script:
        "scripts/postprocess_epinano_differr.py"

# Guitar plot visualization rules
rule guitar_plots:
    input:
        gtf=f"{REFERENCE_DIR}/genes.gtf",
        bed_files=expand(f"{RESULTS_DIR}/{{tool}}/{{sample}}/{{sample}}_{{tool}}_processed.txt", 
                        tool=["CHEUI", "m6Anet", "Nanocompore", "ELIGOS2"], 
                        sample=SAMPLES)
    output:
        mrna_plot=f"{RESULTS_DIR}/summary/guitar_plots_mrna.png",
        ncrna_plot=f"{RESULTS_DIR}/summary/guitar_plots_ncrna.png"
    params:
        sample_name="all_samples"
    conda:
        "envs/r.yaml"
    shell:
        """
        # Create temporary directory with BED files
        temp_bed_dir={RESULTS_DIR}/temp/guitar_bed
        mkdir -p $temp_bed_dir
        
        # Convert processed files to BED format
        for file in {input.bed_files}; do
            if [ -f "$file" ]; then
                # Extract tool and sample names
                tool=$(echo $file | cut -d'/' -f2)
                sample=$(echo $file | cut -d'/' -f3)
                
                # Create BED format (first 6 columns)
                cut -f1-6 "$file" > "$temp_bed_dir/${sample}_${tool}.bed"
            fi
        done
        
        # Generate Guitar plots
        Rscript {workflow.basedir}/scripts/create_guitar_plots.R {input.gtf} $temp_bed_dir {params.sample_name}
        
        # Move plots to final location
        mv {RESULTS_DIR}/guitar_plots/{params.sample_name}/guitar_plot_mrna.png {output.mrna_plot}
        mv {RESULTS_DIR}/guitar_plots/{params.sample_name}/guitar_plot_ncrna.png {output.ncrna_plot}
        
        # Clean up
        rm -rf $temp_bed_dir
        """

# R2Dtool liftover rules
rule r2d_liftover_all:
    input:
        processed_files=expand(f"{RESULTS_DIR}/{{tool}}/{{sample}}/{{sample}}_{{tool}}_processed.txt", 
                              tool=["CHEUI", "DENA", "DRUMMER", "ELIGOS2", "m6Anet", "Nanocompore", "MINES"], 
                              sample=SAMPLES),
        gtf=f"{REFERENCE_DIR}/genes.gtf"
    output:
        summary=f"{RESULTS_DIR}/summary/liftover_summary.tsv",
        done=touch(f"{RESULTS_DIR}/summary/liftover_complete.txt")
    params:
        r2d_tool=config["utilities"]["r2d_tool"]
    shell:
        """
        # Create liftover directory
        liftover_dir={RESULTS_DIR}/liftover
        mkdir -p $liftover_dir
        
        # Process each tool's results
        for file in {input.processed_files}; do
            if [ -f "$file" ]; then
                # Extract tool and sample names
                tool=$(echo $file | cut -d'/' -f2)
                sample=$(echo $file | cut -d'/' -f3)
                
                # Run liftover
                output_file="$liftover_dir/${sample}_${tool}_liftover.txt"
                {params.r2d_tool} liftover -H -g {input.gtf} -i "$file" > "$output_file"
                
                echo "Liftover completed: ${sample}_${tool}" >> {output.summary}
            fi
        done
        
        # Create summary statistics
        echo "Tool\tSample\tLifted_Count\tOriginal_Count" > {output.summary}
        
        for lifted_file in $liftover_dir/*_liftover.txt; do
            if [ -f "$lifted_file" ]; then
                original_file=$(echo $lifted_file | sed 's|liftover/|/|g' | sed 's|_liftover.txt$|_processed.txt|g')
                
                if [ -f "$original_file" ]; then
                    original_count=$(wc -l < "$original_file")
                    lifted_count=$(wc -l < "$lifted_file")
                    
                    # Extract tool and sample
                    filename=$(basename "$lifted_file" _liftover.txt)
                    sample=$(echo $filename | cut -d'_' -f1)
                    tool=$(echo $filename | cut -d'_' -f2-)
                    
                    echo "$tool\t$sample\t$lifted_count\t$original_count" >> {output.summary}
                fi
            fi
        done
        """
