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
        log = "logs/fastp/pe/{BioSample}/{library_type}/{rep}/{Run}.log"
    input:
        fq1 = lambda wildcards: 'raw/' + runTable.loc[(runTable.Run == wildcards['Run']) & (runTable.BioSample == wildcards['BioSample']) & (runTable.replicate == wildcards['rep']) & (runTable.library_type == wildcards['library_type']), ['fq1']].values.tolist()[0][0],
        fq2 = lambda wildcards: 'raw/' + runTable.loc[(runTable.Run == wildcards['Run']) & (runTable.BioSample == wildcards['BioSample']) & (runTable.replicate == wildcards['rep']) & (runTable.library_type == wildcards['library_type']), ['fq2']].values.tolist()[0][0]
    output:
        out1 = "fastp/trimmed/pe/{BioSample}/{library_type}/{rep}/{Run}_1.fastq.gz",
        out2 = "fastp/trimmed/pe/{BioSample}/{library_type}/{rep}/{Run}_2.fastq.gz",
        report_html = "fastp/report/pe/{BioSample}/{library_type}/{rep}/{Run}.fastp.html",
        report_json = "fastp/report/pe/{BioSample}/{library_type}/{rep}/{Run}.fastp.json"
    shell:
        """
            fastp -i {input.fq1} -I {input.fq2}\
                  -o {output.out1} -O {output.out2}\
                  --html {output.report_html} --json {output.report_json}\
                  --length_required 30\
                  --detect_adapter_for_pe\
                  --thread {threads} 2>>{log.log}
        """
