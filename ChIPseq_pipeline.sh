nohup snakemake -j 999 --cluster-config cluster.json \
--cluster "qsub -M {cluster.email} -m {cluster.EmailNotice} -N {cluster.jname} \
-l select=1:ncpus={cluster.cpu} -o {cluster.output} -e {cluster.error}" &> snakemake.log&

