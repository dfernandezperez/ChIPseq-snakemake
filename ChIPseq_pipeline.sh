# Usage --> ./ChIPseq_pipeline.sh       (the samples.json file must be created first, also remember to check the cluster and config files to put everything to your needs
mkdir clusterLogs # Create folder to store the cluster output and error files
nohup snakemake -j 999 --cluster-config cluster.json --latency-wait 120 \
--cluster "qsub -M {cluster.email} -m {cluster.EmailNotice} -N {cluster.jname} \
-l select=1:ncpus={cluster.cpu} -o {cluster.output} -e {cluster.error}" &> snakemake.log&
