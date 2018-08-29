#!/bin/sh
# properties = {"input": ["fastq/E14-K27me3.fastq"], "wildcards": {"sample": "E14-K27me3"}, "params": {"bowtie": "--chunkmbs 1024 -m 1 --best /hpcnfs/techunits/bioinformatics/db/bowtie/mm9/mm9"}, "jobid": 18, "cluster": {"cpu": 10, "email": "ieo4462@ieo.it", "EmailNotice": "a", "jname": "ChIPseq.align.sample=E14-K27me3", "output": "clusterLogs/align.sample=E14-K27me3.out", "error": "clusterLogs/align.sample=E14-K27me3.err"}, "output": ["02aln/E14-K27me3.bam"], "rule": "align", "resources": {}, "local": false, "log": ["00log/E14-K27me3.align", "00log/E14-K27me3.markdup"], "threads": 10, "type": "single"}
cd /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake && \
/hpcnfs/home/ieo4462/.conda/envs/snakemake-tutorial/bin/python \
-m snakemake 02aln/E14-K27me3.bam --snakefile /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.hf65yspr fastq/E14-K27me3.fastq --latency-wait 120 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --mode 2  --allowed-rules align  && touch "/hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.hf65yspr/18.jobfinished" || (touch "/hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.hf65yspr/18.jobfailed"; exit 1)

