#!/bin/sh
# properties = {"input": ["fastq/E14-Ring1b.fastq"], "wildcards": {"sample": "E14-Ring1b"}, "params": {"bowtie": "--chunkmbs 1024 -m 1 --best /hpcnfs/techunits/bioinformatics/db/bowtie/mm9/mm9"}, "jobid": 19, "cluster": {"cpu": 10, "email": "ieo4462@ieo.it", "EmailNotice": "a", "jname": "ChIPseq.align.sample=E14-Ring1b", "output": "clusterLogs/align.sample=E14-Ring1b.out", "error": "clusterLogs/align.sample=E14-Ring1b.err"}, "output": ["02aln/E14-Ring1b.bam"], "rule": "align", "resources": {}, "local": false, "log": ["00log/E14-Ring1b.align", "00log/E14-Ring1b.markdup"], "threads": 10, "type": "single"}
cd /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake && \
/hpcnfs/home/ieo4462/.conda/envs/snakemake-tutorial/bin/python \
-m snakemake 02aln/E14-Ring1b.bam --snakefile /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.hf65yspr fastq/E14-Ring1b.fastq --latency-wait 120 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --mode 2  --allowed-rules align  && touch "/hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.hf65yspr/19.jobfinished" || (touch "/hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.hf65yspr/19.jobfailed"; exit 1)

