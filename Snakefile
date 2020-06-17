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
config['project'] = 'mcf10a_chip'

##### load additional functions #####

include: "scripts/helper.py"

##### global variables/constraints #####
#wildcard_constraints:
#    run="[^_]*",
#    resolution="\d+"

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

rule all_plotCorrelation:
    input:
        expand('deeptools/plotCorrelation/{aggregate}.npz', aggregate = list(runTable.aggregate_column.unique()))


##### load additional workflow rules #####
include: "rules/fastp.smk"
include: "rules/align.smk"
include: "rules/deeptools.smk"
