# Usage --> ./execute_pipeline.sh  <rules> (can be nothing or a rule defined in the snakefile without wildcards in the input files)
mkdir -p .clusterLogs # Create folder to store the cluster output and error files
snakemake -j 1 --cleanup-shadow 
snakemake -j 1 --unlock
nohup snakemake --profile workflow/snakemake_profile "$@" &>> snakemake.log&
