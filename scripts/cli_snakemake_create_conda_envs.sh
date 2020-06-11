/home/150/sxk150/miniconda3/envs/snakemake/bin/snakemake -s /home/150/sxk150/mcf10a_chip/Snakefile ${1}\
    --configfile /home/150/sxk150/mcf10a_chip/config.yaml\
	--use-conda\
	-d /home/150/sxk150/data/mcf10a-hic/ChIP\
	--rerun-incomplete \
        --local-cores 1\
	--cluster-config /home/150/sxk150/mcf10a_chip/cluster.json\
        --keep-going\
	--config machine='gadi'\
	-pr --create-envs-only
