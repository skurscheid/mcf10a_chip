# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-04"


"""
Rules for aligning reads with bowtie2
(https://github.com/BenLangmead/bowtie2)

For usage, include this in your workflow.
"""

def get_index(machine, config):
    """ returns path to index"""
    return config["params"]["bowtie2"]["index"][machine]

singularity: "docker://skurscheid/snakemake_baseimage:0.2"

rule bowtie2_pe_global:
    """ Runs alignment of paired-end fastq files"""
    conda:
        "../envs/alignment.yaml"
    threads:
        12
    params:
        fastq_suffix = ['.end1.fastq.gz', '.end2.fastq.gz'],
        index = get_index(machine, config),
        cli_params_global = config['params']['bowtie2']['cli_params_global'],
        samtools_params_global = "-F 4 -bS"
    log:
        logfile = "logs/bowtie2_global/pe/{BioSample}/{library_type}/{rep}/{Run}.log"
    input:
        rules.run_fastp_pe.output.out1, rules.run_fastp_pe.output.out2
    output:
        bam = temp("bowtie2/align_global/pe/{BioSample}/{library_type}/{rep}/{Run}.bam")
    shell:
        """
            export cli_threads=$(expr {threads} - 2);\
            bowtie2\
                    -x {params.index}\
                    -p $cli_threads\
                    -1 {input[0]} -2 {input[1]}\
                    {params.cli_params_global}\
                    --rg-id BMG\
                    --rg SM:{wildcards.Run}:{wildcards.BioSample}:{wildcards.library_type}:{wildcards.rep}\
                    2>> {log.logfile}\
            | samtools view {params.samtools_params_global} - > {output.bam}
        """
        
rule bam_quality_filter:
    conda:
        "../envs/alignment.yaml"
    version:
        "1.0"
    group:
        "alignment_post"
    log:
        logfile = "logs/samtools/quality_filtered/pe/{BioSample}/{library_type}/{rep}/{Run}.log"
    params:
        qual = config["params"]["general"]["alignment_quality"]
    input:
        rules.bowtie2_pe_global.output
    output:
        temp("samtools/quality_filtered/pe/{BioSample}/{library_type}/{rep}/{Run}.bam")
    shell:
        "samtools view -b -h -q {params.qual} {input} > {output} 2>{log.logfile}"

rule bam_sort:
    conda:
        "../envs/alignment.yaml"
    version:
        "1.0"
    threads:
        4
    group:
        "alignment_post"
    log:
        logfile = "logs/samtools/sort/pe/{BioSample}/{library_type}/{rep}/{Run}.log"
    params:
        temp_dir = 'temp/'
    input:
        rules.bam_quality_filter.output
    output:
        temp("samtools/sort/pe/{BioSample}/{library_type}/{rep}/{Run}.bam")
    shell:
        "samtools sort -@ {threads} {input} -T {params.temp_dir}{wildcards.BioSample}{wildcards.rep}{wildcards.Run}.sorted -o {output}"

rule bam_mark_duplicates:
    conda:
        "../envs/alignment.yaml"
    version:
        "1.0"
    group:
        "alignment_post"
    log:
        logfile = "logs/picardTools/MarkDuplicates/pe/{BioSample}/{library_type}/{rep}/{Run}.log"
    threads:
        4
    params:
        temp = "temp"
    input:
        rules.bam_sort.output
    output:
        out= temp("picardTools/MarkDuplicates/pe/{BioSample}/{library_type}/{rep}/{Run}.bam"),
        metrics = "picardTools/MarkDuplicates/pe/{BioSample}/{library_type}/{rep}/{Run}.metrics.txt"
    shell:
        """
            picard MarkDuplicates -XX:ParallelGCThreads={threads} -Xms2g -Xmx8g\
            INPUT={input}\
            OUTPUT={output.out}\
            ASSUME_SORTED=TRUE\
            METRICS_FILE={output.metrics} 2>{log.logfile}
        """

rule bam_rmdup:
    conda:
        "../envs/alignment.yaml"
    group:
        "alignment_post"
    log:
        logfile = "logs/samtools/rmdup/pe/{BioSample}/{library_type}/{rep}/{Run}.log"
    input:
        rules.bam_mark_duplicates.output.out
    output:
        "samtools/rmdup/pe/{BioSample}/{library_type}/{rep}/{Run}.bam"
    shell:
        "samtools rmdup -s {input} {output} 2>{log.logfile}"

rule bam_index:
    conda:
        "../envs/alignment.yaml"
    group:
        "alignment_post"
    log:
        logfile = "logs/samtools/index/pe/{BioSample}/{library_type}/{rep}/{Run}.log"
    input:
        rules.bam_rmdup.output
    output:
        "samtools/rmdup/pe/{BioSample}/{library_type}/{rep}/{Run}.bam.bai"
    shell:
        "samtools index {input} {output} 2>{log.logfile}"

