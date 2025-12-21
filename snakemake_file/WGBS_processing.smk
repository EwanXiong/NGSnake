'''
Title: Snakefile for Whole Genome Bisulfite Sequencing preprocess

Author: Yifan Xiong

Usage:
    snakemake -s WGBS_preprocessin.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##
import re
workdir: config["workdir"]
bismark_index = config["bismark_index"] # Updated from star_index
samples = config["samples"]

rule all:
    input:
        multiqc = "multiqc/multiqc_report.html",
        bam = expand("bismark/{sample}/{sample}_1_bismark_bt2_pe.bam", sample = samples),
        methy = expand("bismark/{sample}/{sample}.bismark.cov.gz", sample = samples)

rule qc:
    wildcard_constraints: sample = "|".join([re.escape(x) for x in samples])
    container: "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    input:
        R1 = "fastq/{sample}_1.fastq.gz",
        R2 = "fastq/{sample}_2.fastq.gz"
    output:
        R1 = "fastp/{sample}_1.fastq.gz",
        R2 = "fastp/{sample}_2.fastq.gz",
        json = "fastp/log/{sample}.json"
    params: html = "fastp/log/{sample}.html"
    threads: 8
    shell:
        "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} "
        "-w {threads} -h {params.html} -j {output.json}"

rule bismark_mapping:
    wildcard_constraints: sample = "|".join([re.escape(x) for x in samples])
    container: "docker://quay.io/biocontainers/bismark:0.24.2--hdfd78af_0"
    input:
        R1 = "fastp/{sample}_1.fastq.gz",
        R2 = "fastp/{sample}_2.fastq.gz"
    output:
        bam = "bismark/{sample}/{sample}_1_bismark_bt2_pe.bam",
        report = "bismark/{sample}/{sample}_1_bismark_bt2_PE_report.txt"
    threads: 16
    params: outdir = "bismark/{sample}/"
    shell:
        "bismark --bowtie2 -p {threads} --output_dir {params.outdir} "
        "{bismark_index} -1 {input.R1} -2 {input.R2}"

rule bismark_methylation_extract:
    container: "docker://quay.io/biocontainers/bismark:0.24.2--hdfd78af_0"
    input: "bismark/{sample}/{sample}_1_bismark_bt2_pe.bam"
    output:
        cov = "bismark/{sample}/{sample}.bismark.cov.gz"
    params: outdir = "bismark/{sample}/"
    threads: 4
    shell:
        "bismark_methylation_extractor --paired-end --bedGraph --counts "
        "--gzip -o {params.outdir} {input}"

rule multiqc:
    container: "docker://multiqc/multiqc:v1.28"
    input:
        json = expand("fastp/log/{sample}.json", sample = samples),
        reports = expand("bismark/{sample}/{sample}_1_bismark_bt2_PE_report.txt", sample = samples)
    output: "multiqc/multiqc_report.html"
    shell: "multiqc . -o multiqc/"

## DON'T CHANGE ABOVE THIS LINE ##
