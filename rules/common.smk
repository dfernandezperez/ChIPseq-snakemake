def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])


def is_spike(sample):
    return SAMPLES.SPIKE[sample] == True


def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("fastq/{sample}.{group}.fastq", group=[1, 2], **wildcards)
    # single end sample
    return "fastq/{sample}.se.fastq".format(**wildcards)


def get_trimmed_forward(wildcards):
	if not is_single_end(**wildcards):
	    # paired-end sample
	    return "fastq/{sample}.1.fastq".format(**wildcards)
	# single end sample
	return "fastq/{sample}.se.fastq".format(**wildcards)


def get_trimmed_spike(wildcards):
    if is_spike(**wildcards):
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("fastq/{sample}.{group}.fastq", group=[1, 2], **wildcards)
        # single end sample
        return "fastq/{sample}.se.fastq".format(**wildcards)


def get_bam(wildcards):
    if not is_spike(wildcards.sample):
        return "02aln/{sample}.bam".format(wildcards.sample)
    return "02aln/{sample}.clean.bam".format(wildcards.sample)


def get_bam_cntrl(wildcards):
    if not is_spike(wildcards.sample):
        return { "case": "02aln/{sample}.clean.bam".format(sample=wildcards.sample),
            "reference": "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/{genome}/input_{control}.bam".format(
            genome = SAMPLES.GENOME[wildcards.sample], 
            control = wildcards.control) }

    return { "case": "02aln/{sample}.clean.bam".format(sample=wildcards.sample),
             "reference": "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/{genome}/input_{control}.bam".format(
             genome = SAMPLES.GENOME[wildcards.sample], 
             control = wildcards.control),
             "spike": "02aln_dm/{sample}_spike.clean.bam".format(sample=wildcards.sample) }

