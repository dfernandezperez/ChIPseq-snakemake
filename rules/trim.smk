def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()

rule cp_fastq_pe:
    input:
        get_fastq
    output:
        fastq1=temp("/hpcnfs/scratch/temporary/DP/fastq/{sample}-{lane}.1.fastq.gz"),
        fastq2=temp("/hpcnfs/scratch/temporary/DP/fastq/{sample}-{lane}.2.fastq.gz")
    message:
        "Copying fastq files {input}"
    shell:
        """
        ln -s {input[0]} {output.fastq1}
        ln -s {input[1]} {output.fastq2}
        """


rule cp_fastq_se:
    input:
        get_fastq
    output:
        temp("/hpcnfs/scratch/temporary/DP/fastq/{sample}-{lane}.fastq.gz"),
    message:
        "Copying fastq files {input}"
    shell:
        """
        ln -s {input} {output}
        """


rule fastp_pe:
	input:
		fw = lambda w: expand("/hpcnfs/scratch/temporary/DP/fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=units.loc[w.sample].itertuples()),
		rv = lambda w: expand("/hpcnfs/scratch/temporary/DP/fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=units.loc[w.sample].itertuples())
	output:
		fastq1 = temp("/hpcnfs/scratch/temporary/DP/fastq/{sample}.1.fastq"),
		fastq2 = temp("/hpcnfs/scratch/temporary/DP/fastq/{sample}.2.fastq")
	log:
		"00log/fastp/{sample}.log"
	threads:
		CLUSTER["fastp_pe"]["cpu"]
	params:
		fastp_params = config["params"]["fastp"]["pe"],
		tmp_fw       = "/hpcnfs/scratch/temporary/DP/fastq/{sample}.1.fastq.tmp.gz",
		tmp_rv       = "/hpcnfs/scratch/temporary/DP/fastq/{sample}.2.fastq.tmp.gz"
	message:
		"Processing fastq files from {input}"
	shadow:
		"minimal"
	benchmark:
		".benchmarks/{sample}.merge_fastqs.benchmark.txt"
	shell:
		"""
		cat {input.fw} > {params.tmp_fw}
		cat {input.rv} > {params.tmp_rv}
		fastp -i {params.tmp_fw} \
		-I {params.tmp_rv} \
		-o {output.fastq1} \
		-O {output.fastq2} \
		-w {threads} \
		{params.fastp_params} 2> {log}
		"""


rule fastp_se:
	input:
		lambda w: expand("/hpcnfs/scratch/temporary/DP/fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=units.loc[w.sample].itertuples()),
	output:
		temp("/hpcnfs/scratch/temporary/DP/fastq/{sample}.se.fastq")
	log:
		"00log/fastp/{sample}.log"
	threads:
		CLUSTER["fastp_se"]["cpu"]
	params:
		fastp_params = config["params"]["fastp"]["se"],
	message:
		"Processing fastq files from {input}"
	shadow:
		"minimal"
	benchmark:
		".benchmarks/{sample}.merge_fastqs.benchmark.txt"
	shell:
		"""
		zcat {input} | \
		fastp -o {output} \
		-w {threads} {params.fastp_params} 2> {log}
		"""
