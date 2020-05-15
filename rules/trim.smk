def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()

rule cp_fastq_pe:
    input:
        get_fastq
    output:
        fastq1=temp("fastq/{sample}-{lane}.1.fastq.gz"),
        fastq2=temp("fastq/{sample}-{lane}.2.fastq.gz")
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
        temp("fastq/{sample}-{lane}.fastq.gz"),
    message:
        "Copying fastq files {input}"
    shell:
        """
        ln -s {input} {output}
        """


rule mergeFastq_pe:
	input:
		fw = lambda w: expand("fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=units.loc[w.sample].itertuples()),
		rv = lambda w: expand("fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=units.loc[w.sample].itertuples())
	output:
		fastq1 = temp("fastq/{sample}.1.fastq"),
		fastq2 = temp("fastq/{sample}.2.fastq")
	log:
		"00log/fastp/{sample}.log"
	message:
		"Processing fastq files from {input}"
	shell:
		"""
		cat {input.fw} > {output.fastq1}
		cat {input.rv} > {output.fastq1}
		"""


rule mergeFastq_se:
	input:
		lambda w: expand("fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=units.loc[w.sample].itertuples()),
	output:
		temp("fastq/{sample}.se.fastq")
	log:
		"00log/fastp/{sample}.log"
	message:
		"Processing fastq files from {input}"
	shell:
		"""
		cat {input} > {output}
		"""