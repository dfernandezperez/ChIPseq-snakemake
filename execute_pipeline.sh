# Usage --> ./execute_pipeline.sh  <rules> (can be nothing or a rule defined in the snakefile without wildcards in the input files)
snakemake --unlock
mkdir -p .clusterLogs # Create folder to store the cluster output and error files
nohup snakemake --profile ./snakemake_profile "$@" &>> snakemake.log&