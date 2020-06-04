# Usage --> ./execute_pipeline.sh  <rules> (can be nothing or a rule defined in the snakefile without wildcards in the input files)
mkdir -p results/.clusterLogs
snakemake -j 1 --unlock
nohup snakemake --profile workflow/snakemake_profile "$@" &>> results/snakemake.log&
