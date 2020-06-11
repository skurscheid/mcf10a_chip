from pathlib import Path
import os
import yaml
import pandas as pd

# functions used in Snakefile
def make_targets_from_runTable(runTable):
    t = []
    for index, row in runTable.iterrows():
        e = list(row[['BioSample', 'library_type', 'replicate', 'Run']])
        p = "/".join(e)
        t.append(p)
    return(t)

def create_testing_input(base_path, units):
    """creates test files for snakemake run"""
    for index, row in units.iterrows():
        fq1, fq2 = Path(base_path ,row['fq1']), Path(base_path ,row['fq2'])
        p = Path(os.path.split(fq1)[0])
        p.mkdir(parents = True, exist_ok = True)
        fq1.touch(exist_ok = True)
        fq2.touch(exist_ok = True)

# functions for fastp rules
def fastp_input(wildcards, runTable):
    """function for creating gathering input files for fastp processing"""
    return(runTable.loc[(runTable.Run == wildcards['Run']) & (runTable.BioSample == wildcards['BioSample']) & (runTable.replicate == wildcards['rep']) & (runTable.library_type == wildcards['library_type']), ['fq1', 'fq2']].values.tolist())
