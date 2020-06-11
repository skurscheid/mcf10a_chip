# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

import pandas as pd
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####

min_version("5.1.2")

##### load codMCF10AshZsheets #####

configfile: "config.yaml"

runTable = pd.read_csv("ChIPRunTable.csv", sep = ",")
machine = config['machine']

##### load additional functions #####

include: "scripts/helper.py"

##### load additional workflow rules #####
include: "rules/fastp.smk"
include: "rules/align.smk"

##### global variables/constraints #####
wildcard_constraints:
    run="[^_]*",
    resolution="\d+"

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
        unpack(fastp_input)
        #fq1 = 'raw/SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq/H2AZ_TGFb_10A_high_rep2/Florida_10_S5_R1_001.fastq.gz',
        #fq2 = 'raw/SN501_0087_DTremethick_JCSMR_MCF10A_ChIPSeq/H2AZ_TGFb_10A_high_rep2/Florida_10_S5_R2_001.fastq.gz'
    output:
        out1 = "fastp/trimmed/pe/{biosample}/{library_type}/{replicate}/{run}_1.fastq.gz",
        out2 = "fastp/trimmed/pe/{biosample}/{library_type}/{replicate}/{run}_2.fastq.gz",
        report_html = "fastp/report/pe/{biosample}/{library_type}/{replicate}/{run}.fastp.html",
        report_json = "fastp/report/pe/{biosample}/{library_type}/{replicate}/{run}.fastp.json"
    shell:
        """
            fastp -i {input.fq1} -I {input.fq2}\
                  -o {output.out1} -O {output.out2}\
                  --html {output.report_html} --json {output.report_json}\
                  --length_required 30\
                  --detect_adapter_for_pe\
                  --thread {threads} 2>>{log.log}
        """


##### build targets #####
#rule all:
#    input:
        # The first rule should define the default target files
        # Subsequent target rules can be specified below. They should start with all_*.

rule all_trim:
    input:
        expand("fastp/trimmed/pe/{file}{end}.fastq.gz",
                file = make_targets_from_runTable(runTable),
                end = [config["params"]["general"]["end1_suffix"], config["params"]["general"]["end2_suffix"]])

rule all_align:
    input:
        expand("samtools/rmdup/pe/{file}.bam.bai", file = make_targets_from_runTable(runTable))


