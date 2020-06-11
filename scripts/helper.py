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
def fastp_input(runTable, wildcards):
    """function for creating gathering input files for fastp processing"""
    t = []
    row = runTable.loc[(runTable.Run == wildcards['Run']) & (runTable.BioSample == wildcards['BioSample']) & (runTable.replicate == wildcards['rep']) & (runTable.library_type == wildcards['library_type'])]
    fq1, fq2 = row['fq1'].to_string(index=False), row['fq2'].to_string(index=False)
    t.append(fq1.strip())
    t.append(fq2.strip())
    return(t)

def fastp_targets(units):
    """function for creating snakemake targets for executing fastp rule"""
    t = []
    for index, row in units.iterrows():
        t.append(row['batch'] + "/" + row['sample_id'] + "_" + row['lane'] + "_" + str(row['replicate']))
    return(t)
