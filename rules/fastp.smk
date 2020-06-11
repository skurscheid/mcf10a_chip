__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for trimming reads with fastq
(https://github.com/OpenGene/fastp)

For usage, include this in your workflow.
"""

rule run_fastp_pe:
    conda:
        "../envs/fastqProcessing.yaml"
    version:
        2
    threads:
        4
    log:
        log = "logs/fastp/pe/{biosample}/{library_type}/{replicate}/{run}.log"
    input:
        fastp_input
    output:
        out1 = "fastp/trimmed/pe/{biosample}/{library_type}/{replicate}/{run}_1.fastq.gz",
        out2 = "fastp/trimmed/pe/{biosample}/{library_type}/{replicate}/{run}_2.fastq.gz",
        report_html = "fastp/report/pe/{biosample}/{library_type}/{replicate}/{run}.fastp.html",
        report_json = "fastp/report/pe/{biosample}/{library_type}/{replicate}/{run}.fastp.json"
    shell:
        """
            fastp -i {input[0]} -I {input[1]}\
                  -o {output.out1} -O {output.out2}\
                  --html {output.report_html} --json {output.report_json}\
                  --length_required 30\
                  --detect_adapter_for_pe\
                  --thread {threads} 2>>{log.log}
        """
