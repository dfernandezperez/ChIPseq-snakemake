#!/bin/sh
# properties = {"params": {"bowtie_mm": "--chunkmbs 1024 -m 1 --best /hpcnfs/techunits/bioinformatics/db/bowtie/mm9/mm9", "bowtie_dm": "--chunkmbs 1024 -m 1 --best /hpcnfs/data/DP/dm6/dm6"}, "wildcards": {"sample": "E14-K27me3"}, "input": ["fastq/E14-K27me3.fastq"], "jobid": 18, "rule": "align_spike", "cluster": {"cpu": 1, "email": "ieo4462@ieo.it", "EmailNotice": "a", "jname": "ChIPseq.align_spike.sample=E14-K27me3", "output": "clusterLogs/align_spike.sample=E14-K27me3.out", "error": "clusterLogs/align_spike.sample=E14-K27me3.err"}, "resources": {}, "type": "single", "output": ["02aln/E14-K27me3.bam", "02aln_dm6/E14-K27me3_dm6.bam"], "local": false, "log": ["00log/E14-K27me3.align", "00log/E14-K27me3.dm6_align", "00log/E14-K27me3.markdup"], "threads": 10}
cd /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake && \
/hpcnfs/home/ieo4462/.conda/envs/snakemake-tutorial/bin/python \
-m snakemake 02aln/E14-K27me3.bam --snakefile /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.owjaxogz fastq/E14-K27me3.fastq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://bitbucket.org/snakemake/snakemake-wrappers/raw/ \
   --nocolor \
--notemp --no-hooks --nolock --mode 2  --allowed-rules align_spike  && touch "/hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.owjaxogz/18.jobfinished" || (touch "/hpcnfs/scratch/DP/dfernand/Test3/ChIPseq-snakemake/.snakemake/tmp.owjaxogz/18.jobfailed"; exit 1)

