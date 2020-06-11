#!/bin/bash
#PBS -P zi26
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -q express
#PBS -e /home/150/sxk150/qsub_error
#PBS -o /home/150/sxk150/qsub_out
#PBS -l ncpus=4
#PBS -l mem=32GB
#PBS -M skurscheid@gmail.com
#PBS -m abe
#PBS -l storage=scratch/kv78+gdata/kv78

target=${cli_target}

source ~/.bashrc

/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_chip/Snakefile ${cli_target}\
    --configfile /home/150/sxk150/mcf10a_chip/config.yaml\
	--use-conda\
	--cluster "qsub -P {cluster.P}\
                    -l ncpus={threads} \
                    -q {cluster.queue} \
                    -l mem={cluster.mem} \
                    -l wd\
                    -l walltime={cluster.walltime}\
		            -l storage={cluster.storage}\
		            -l jobfs={cluster.jobfs}\
                    -M {cluster.M}\
                    -m {cluster.m}\
                    -e {cluster.error_out_dir} \
                    -o {cluster.std1_out_dir}" \
	--jobs 100\
	-d /home/150/sxk150/data/mcf10a-hic/ChIP\
	--rerun-incomplete \
    --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10a_chip/cluster.json\
    --config machine=gadi\
    --keep-going\
	-pr 

