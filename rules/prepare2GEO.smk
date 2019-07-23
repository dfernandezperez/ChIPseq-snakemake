rule geo_fastp_pe:
	input:
		fw = lambda w: expand("fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=units.loc[w.sample].itertuples()),
		rv = lambda w: expand("fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=units.loc[w.sample].itertuples())
	output:
		fastq1 = "GEO/fastq/{sample}.1.fastq.gz",
		fastq2 = "GEO/fastq/{sample}.2.fastq.gz"
	log:
		"00log/GEO_fastp/{sample}.log"
	threads:
		CLUSTER["fastp_pe"]["cpu"]
	params:
		fastp_params = config["params"]["fastp"]["pe"],
		tmp_fw       = "GEO/fastq/{sample}.1.fastq.tmp.gz",
		tmp_rv       = "GEO/fastq/{sample}.2.fastq.tmp.gz"
	message:
		"Processing fastq files from {input}"
	shadow:
		"minimal"
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


rule geo_fastp_se:
	input:
		lambda w: expand("fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=units.loc[w.sample].itertuples()),
	output:
		"GEO/fastq/{sample}.se.fastq.gz"
	log:
		"00log/GEO_fastp/{sample}.log"
	threads:
		CLUSTER["fastp_se"]["cpu"]
	params:
		fastp_params = config["params"]["fastp"]["se"],
	message:
		"Processing fastq files from {input}"
	shadow:
		"minimal"
	shell:
		"""
		fastp -i {input} \
        -o {output} \
		-w {threads} \
        {params.fastp_params} 2> {log}
		"""

rule geo_peaks:
	input:
		"03peak_macs2/{sample}_{control}-input/{sample}_peaks_p{pvalue}.bed"
	output:
		"GEO/peaks/{sample}_{control}-input_peaks_p{pvalue}.bed"
	log:
		"00log/GEO/peaks/{sample}_{control}_{pvalue}.log"
	shell:
		"""
		cp {input} {output} 2> {log}
		"""

#----------------- MD5sum generation ----------------#		
rule md5sum_peaks:
	input:
		expand(
			"GEO/peaks/{sample}_{control}-input_peaks_p" + config["params"]["macs2"]["filt_peaks_pval"] + ".bed",
			zip,
			sample  = ALL_SAMPLES,
			control = ALL_CONTROLS
		)
	output:
		"GEO/md5sum/md5sum_peaks.txt"
	log:
		"00log/GEO/peaks/md5sum.log"
	shell:
		"""
		md5sum {input} > {output} 2> {log}
		"""

rule md5sum_fastqs:
	input:
		expand("GEO/fastq/{sample}.se.fastq.gz", sample = ALL_SAMPLES_SE),
		expand(["GEO/fastq/{sample}.1.fastq.gz", "GEO/fastq/{sample}.2.fastq.gz"], sample = ALL_SAMPLES_PE)
	output:
		"GEO/md5sum/md5sum_fastqs.txt"
	log:
		"00log/GEO/fastqs/md5sum.log"
	shell:
		"""
		md5sum {input} > {output} 2> {log}
		"""