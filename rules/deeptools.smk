__author__ = "Sebastian Kurscheid (sebastian.kurscheid@anu.edu.au)"
__license__ = "MIT"
__date__ = "2020-02-04"

# vim: syntax=python tabstop=4 expandtab
# coding: utf-8


"""
Rules for processing HTS data with deepTools
(https://deeptools.readthedocs.io/en/develop/)

For usage, include this in your workflow.
"""

# input & other functions
def get_multi_bam_summary_input(wildcards):
    aggregate = wildcards['aggregate']
    l = []
    sel_rows = runTable['aggregate_column'] == aggregate
    for index, row in runTable[sel_rows].iterrows():
        l.append('/'.join(["samtools/rmdup/pe", row.BioSample, row.library_type, row.replicate, row.Run]) + '.bam')
    return(l)

def get_multi_bam_summary_labels(wildcards):
    aggregate = wildcards['aggregate']
    l = []
    sel_rows = runTable['aggregate_column'] == aggregate
    for index, row in runTable[sel_rows].iterrows():
        l.append('_'.join([row.aggregate_column, row.Run, row.replicate]))
    return(l)

def get_bigwigCompare_inputs(wildcards):
    selected_columns = config['params']['general']['runTable'][config['project']]['selected_columns']
    library_type = config['library_type']
    cell_line = wildcards['cell_line']
    chip_antibody = wildcards['chip_antibody']
    l = []
    selection = runTable[(runTable['aggregate_column'] == chip_antibody) & (runTable[selected_columns[0]] == cell_line)]
    for index, row in selection.iterrows():
        l.append('/'.join(["deeptools/bamCoverage", cell_line, row.aggregate_column, library_type, row.Run]) + '.bw')
    return(l)

# actual rules
rule macs2_predictd:
    version:
        1
    conda:
        '../envs/macs2.yaml'
    threads:
        1
    group:
       'deeptools'
    log:
        logfile = 'logs/macs2/predictd/{cell_line}/{chip_antibody}/{library_type}/{run}.log'
    params:
        gsize = 'hs'
    input:
        'samtools/rmdup/{cell_line}/{chip_antibody}/{library_type}/{run}.bam'
    output:
        file = ('macs2/predictd/{cell_line}/{chip_antibody}/{library_type}/{run}_predictd.R')
    shell:
        """
            macs2 predictd -i {input}\
                           --gsize {params.gsize}\
                           --rfile {output.file} 2>{log.logfile}
        """

rule deeptools_multiBamSummary:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        8
    group:
        "deeptools"
    params:
        labels = get_multi_bam_summary_labels
    log:
        logfile = "logs/deeptools_multiBamSummary/{aggregate}.log"
    input:
        get_multi_bam_summary_input
    output:
        npz = "deeptools/multiBamSummary/{aggregate}.npz"
    shell:
        """
            multiBamSummary bins --bamfiles {input}\
                                 --numberOfProcessors {threads}\
                                 --labels {params.labels}\
                                 --outFileName {output.npz} 2>{log.logfile}
        """

rule deeptools_plotCorrelation:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        1
    group:
        "deeptools"
    params:
    log:
        logfile = "logs/deeptools_plotCorrelation/{cell_line}.log"
    input:
        rules.deeptools_multiBamSummary.output.npz
    output:
        png = "deeptools/plotCorrelation/{cell_line}_heatmap.png"
    shell:
        """
            plotCorrelation -in {input} --whatToPlot heatmap --corMethod pearson -o {output}
        """

rule deeptools_bamCoverage:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        8
    group:
        "deeptools"
    params:
        extendReads = 200,
        binSize = config['params']['deeptools']['binSize'],
        smoothLength = config['params']['deeptools']['smoothLength'],
        effectiveGenomeSize = config['params']['deeptools']['genome_size']['GRCh37_hg19_UCSC'],
        normalizeUsing = config['params']['deeptools']['normalizeUsing']
    log:
        logfile = "logs/deeptools_bamCoverage/{cell_line}/{chip_antibody}/{library_type}/{run}.log"
    input:
        rules.bam_rmdup.output
    output:
        "deeptools/bamCoverage/{cell_line}/{chip_antibody}/{library_type}/{run}.bw"
    shell:
        """
            bamCoverage -b {input}\
                        --numberOfProcessors {threads}\
                        --effectiveGenomeSize {params.effectiveGenomeSize}\
                        --normalizeUsing {params.normalizeUsing}\
                        --extendReads {params.extendReads}\
                        --binSize {params.binSize}\
                        --smoothLength {params.smoothLength}\
                        --outFileName {output} 2>{log.logfile}
        """

rule merge_bigwigs:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        4
    group:
        "deeptools"
    params:
    log:
        logfile = "logs/merge_bigwigs/{cell_line}/{chip_antibody}_coverage.log"
    input:
        chip = get_bigwigCompare_inputs
    output:
        "deeptools/merge_bigwigs/{cell_line}/{chip_antibody}_coverage.bw"
    script:
        "../scripts/deeptools_merge_bigwigs.py"


rule deeptools_bigwigCompare:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        8
    group:
        "deeptools"
    params:
        binSize = config['params']['deeptools']['binSize'],
    log:
        logfile = "logs/merge_bigwigs/{cell_line}/{chip_antibody}_coverage.log"
    input:
        chip = "deeptools/merge_bigwigs/{cell_line}/{chip_antibody}_coverage.bw", 
        input = "deeptools/merge_bigwigs/{cell_line}/Input_coverage.bw" 
    output:
        "deeptools/bigwigCompare/{cell_line}/{chip_antibody}_vs_Input.bw"
    shell:
        """
            bigwigCompare --bigwig1 {input.chip}\
                          --bigwig2 {input.input}\
                          --operation log2\
                          --binSize {params.binSize}\
                          --numberOfProcessors {threads}\
                          --outFileName {output}
        """

rule deeptools_computeMatrix_referencepoint:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        16
    group:
        "deeptools"
    params:
        referencePoint = config['params']['deeptools']['referencePoint'],
        beforeRegionStartLength = config['params']['deeptools']['beforeRegionStartLength'],
        afterRegionStartLength = config['params']['deeptools']['afterRegionStartLength'],
        sortRegions = config['params']['deeptools']['sortRegions'],
        regionsFileName = lambda wildcards: expand("{dir}{figure}", dir = config['params']['deeptools']['annotation_dir'][machine], figure = config['params']['deeptools']['regionsFiles'][wildcards['figure']]),
        annotationDir = config['params']['deeptools']['annotation_dir'][machine]
    log:
        logfile = "logs/deeptools_computeMatrix/{cell_line}/{chip_antibody}_{figure}_matrix.log"
    input:
        rules.deeptools_bigwigCompare.output
    output:
        "deeptools/computeMatrix_referencepoint/{cell_line}/{chip_antibody}_{figure}_matrix.gz"
    shell:
        """
            computeMatrix reference-point --regionsFileName {params.regionsFileName}\
                                          --scoreFileName {input}\
                                          --outFileName {output}\
                                          --referencePoint {params.referencePoint}\
                                          --beforeRegionStartLength {params.beforeRegionStartLength}\
                                          --afterRegionStartLength {params.afterRegionStartLength}\
                                          --sortRegions {params.sortRegions}\
                                          --numberOfProcessors {threads} 2>{log.logfile}
        """

rule deeptools_plotProfile:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        1
    group:
        "deeptools"
    log:
        logfile = "logs/deeptools_plotProfile/{cell_line}/{figure}/{chip_antibody}.log"
    params:
        numPlotsPerRow = config['params']['deeptools']['plotProfile']['numPlotsPerRow'],
        colors = config['params']['deeptools']['plotProfile']['colors'],
        plotType = 'se',
        regionsLabel = config['params']['deeptools']['regionsLabel'],
        yMin = -2.5,
        yMax = 4.5
    input:
        rules.deeptools_computeMatrix_referencepoint.output
    output:
        "deeptools/plotProfile/{cell_line}/{figure}/{chip_antibody}.pdf"
    shell:
        """
            plotProfile -m {input}\
                        -o {output}\
                        --perGroup\
                        --numPlotsPerRow {params.numPlotsPerRow}\
                        --colors {params.colors}\
                        --plotType {params.plotType}\
                        --regionsLabel {params.regionsLabel}\
                        --plotTitle '{wildcards.figure} - {wildcards.chip_antibody}' 2>{log.logfile}
        """

rule deeptools_plotHeatmap:
    version:
        1
    conda:
        "../envs/deeptools.yaml"
    threads:
        1
    group:
        "deeptools"
    log:
        logfile = "logs/deeptools_plotHeatmap/{cell_line}/{chip_antibody}_{figure}.log"
    params:
        regionsLabel = config['params']['deeptools']['regionsLabel'],
        sortRegions = 'keep'
    input:
        rules.deeptools_computeMatrix_referencepoint.output
    output:
        "deeptools/plotHeatmap/{cell_line}/{figure}/{chip_antibody}.pdf"
    shell:
        '''
            plotHeatmap -m {input}\
                        -o {output}\
                        --regionsLabel {params.regionsLabel}\
                        --sortRegions {params.sortRegions}\
                        --plotTitle {wildcards.chip_antibody} 2>{log.logfile}
        '''

