rule geo_mergeFastq_pe:
	input:
		fw = lambda w: expand("fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=units.loc[w.sample].itertuples()),
		rv = lambda w: expand("fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=units.loc[w.sample].itertuples())
	output:
		fastq1 = "GEO/fastq/{sample}.1.fastq.gz",
		fastq2 = "GEO/fastq/{sample}.2.fastq.gz"
	log:
		"00log/GEO_fastp/{sample}.log"
	message:
		"Merging fastq files from {input}"
	shell:
		"""
		cat {input.fw} > {output.fastq1}
		cat {input.rv} > {output.fastq1}
		"""


rule geo_mergeFastq_se:
	input:
		lambda w: expand("fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=units.loc[w.sample].itertuples()),
	output:
		"GEO/fastq/{sample}.se.fastq.gz"
	log:
		"00log/GEO_fastp/{sample}.log"
	message:
		"Processing fastq files from {input}"
	shell:
		"""
		cat {input} > {output}
		"""

rule geo_peaks:
	input:
		"03peak_macs2/{sample}_{control}/{sample}_peaks_p{pvalue}.bed"
	output:
		"GEO/peaks/{sample}_{control}_peaks_p{pvalue}.bed"
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
			"GEO/peaks/{sample}_{control}_peaks_p" + config["params"]["macs2"]["filt_peaks_pval"] + ".bed",
			zip,
			sample  = ALL_IP,
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
		expand("GEO/fastq/{sample}.se.fastq.gz", sample = ALL_IP_SE),
		expand(["GEO/fastq/{sample}.1.fastq.gz", "GEO/fastq/{sample}.2.fastq.gz"], sample = ALL_IP_PE)
	output:
		"GEO/md5sum/md5sum_fastqs.txt"
	log:
		"00log/GEO/fastqs/md5sum.log"
	shell:
		"""
		md5sum {input} > {output} 2> {log}
		"""