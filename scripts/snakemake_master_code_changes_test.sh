#!/bin/bash
#PBS -P pb97
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -q biodev
#PBS -e /home/150/sxk150/qsub_error
#PBS -o /home/150/sxk150/qsub_out
#PBS -l ncpus=1
#PBS -l mem=16GB
#PBS -M skurscheid@gmail.com
#PBS -m abe
#PBS -l storage=scratch/kv78+gdata/kv78

target=${cli_target}

source ~/.bashrc

/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_chip/Snakefile\
    -R `/home/150/sxk150/mcf10a_chip/scripts/cli_snakemake.sh ${target} --lc`\
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
	-d ~/data/mcf10a-hic/ChIP\
	--rerun-incomplete \
    --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10a_chip/cluster.json\
    --config machine=gadi\
    --keep-going\
	-prn

