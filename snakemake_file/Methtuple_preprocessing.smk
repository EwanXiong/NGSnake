'''
Title: Snakefile for Methtuple preprocessing

Author: Yifan Xiong

Usage:
    snakemake -s Methtuple_preprocessing.smk --configfile <path to config.yaml> \
              --cores <int> --use-singularity
'''

## DON'T CHANGE BELOW THIS LINE ##
import re

workdir: config["workdir"]
raw_data = config["raw_data"]
bismark_index = config["bismark_index"]
samples = config["samples"]
threads_set = config["threads"]

sample_constraint = "|".join([re.escape(x) for x in samples])

rule all:
    input:
        multiqc = "QC/multiqc/multiqc_report.html",
        bam = expand("bismark/{sample}/{sample}_bismark_bt2_query_sorted.bam", sample=samples)

rule trimming:
    wildcard_constraints: 
        sample = sample_constraint
    container: 
        "docker://quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    input:
        R1 = os.path.join(raw_data, "{sample}_1.fastq.gz"),
        R2 = os.path.join(raw_data, "{sample}_2.fastq.gz")
    output:
        R1 = "fastp/{sample}_trimmed_1.fastq.gz",
        R2 = "fastp/{sample}_trimmed_2.fastq.gz",
        json = "fastp/log/{sample}.json",
        html = "fastp/log/{sample}.html"
    threads: threads_set
    log: "logs/fastp/{sample}.log"
    shell:
        "fastp -i {input.R1} -I {input.R2} "
        "-o {output.R1} -O {output.R2} "
        "-w {threads} -h {output.html} -j {output.json} "
        "2> {log}"

rule fastqc_clean:
    wildcard_constraints: 
        sample = sample_constraint
    container: 
        "docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    input:
        R1 = "fastp/{sample}_trimmed_1.fastq.gz",
        R2 = "fastp/{sample}_trimmed_2.fastq.gz"
    output:
        html_R1 = "QC/fastqc_clean/{sample}_trimmed_1_fastqc.html",
        zip_R1 = "QC/fastqc_clean/{sample}_trimmed_1_fastqc.zip",
        html_R2 = "QC/fastqc_clean/{sample}_trimmed_2_fastqc.html",
        zip_R2 = "QC/fastqc_clean/{sample}_trimmed_2_fastqc.zip"
    params:
        outdir = "QC/fastqc_clean"
    threads: 2
    log: "logs/QC/fastqc_clean/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -o {params.outdir} -t {threads} {input.R1} {input.R2} 2> {log}
        """

rule multiqc:
    container: 
        "docker://multiqc/multiqc:v1.28"
    input:
        clean_fastqc = expand("QC/fastqc_clean/{sample}_trimmed_{read}_fastqc.zip", 
                             sample=samples, read=["1", "2"]),
        fastp_json = expand("fastp/log/{sample}.json", sample=samples)
    output: 
        "QC/multiqc/multiqc_report.html"
    params:
        search_dir = "QC/",
    log: 
        "logs/QC/multiqc/multiqc.log"
    shell: 
        """
        multiqc {params.search_dir} -o QC/multiqc -f 2> {log}
        """

rule bismark_mapping:
    wildcard_constraints: 
        sample = sample_constraint
    container: 
        "docker://quay.io/biocontainers/bismark:0.24.2--hdfd78af_0"
    input:
        R1 = "fastp/{sample}_trimmed_1.fastq.gz",
        R2 = "fastp/{sample}_trimmed_2.fastq.gz"
    output:
        sam = temp("bismark/{sample}/{sample}_bismark_bt2_pe.sam")
    threads: 32
    params: 
        outdir = "bismark/{sample}/",
        index = bismark_index
    log: "logs/bismark/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        bismark --bowtie2 -p {threads} \
            --output_dir {params.outdir} \
            --basename {wildcards.sample} \
            {params.index} \
            -1 {input.R1} -2 {input.R2} 2> {log}
        """

rule Sort_by_Queryname:
    container:
        "docker://quay.io/biocontainers/samtools:1.19.2--h50ea8bc_0"
    input:
        sam = "bismark/{sample}/{sample}_bismark_bt2_pe.sam"
    output:
        qsort_bam = "bismark/{sample}/{sample}_bismark_bt2_query_sorted.bam"
    threads: 32
    shell:
        "samtools sort -n -@ {threads} -o {output.qsort_bam} {input.sam}"


## DON'T CHANGE ABOVE THIS LINE ##