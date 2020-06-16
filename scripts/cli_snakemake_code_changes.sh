/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_chip/Snakefile\
    --configfile /home/150/sxk150/mcf10a_chip/config.yaml\
	--use-conda\
	-d /home/150/sxk150/data/mcf10a-hic/ChIP\
	--rerun-incomplete \
        --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10a_chip/cluster.json\
        --keep-going\
	-pr\
	--config machine='gadi'\
	-R `/home/150/sxk150/mcf10a_chip/scripts/cli_snakemake.sh ${1} --lc`\
	${2}
