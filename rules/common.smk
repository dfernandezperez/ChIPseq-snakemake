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


def get_trimmed_spike(wildcards):
    if is_spike(**wildcards):
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("fastq/{sample}.{group}.fastq", group=[1, 2], **wildcards)
        # single end sample
        return "fastq/{sample}.se.fastq".format(**wildcards)


# Used to use just the forward reads to create the fastqc report
def get_trimmed_forward(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return "fastq/{sample}.1.fastq".format(**wildcards)
    # single end sample
    return "fastq/{sample}.se.fastq".format(**wildcards)


def get_bam(wildcards):
    if not is_spike(**wildcards):
        return "02aln/{sample}.bam.tmp".format(**wildcards)
    return "02aln/{sample}.bam.tmp.clean".format(**wildcards)


def set_read_extension(wildcards):
    if is_single_end(wildcards.sample):
        return "--extReads " + str(config['read_extension'])
    return "--extReads"

def get_bam_cntrl(wildcards):
    if not is_spike(wildcards.sample):
        return { "case": "02aln/{sample}.bam".format(sample=wildcards.sample),
            "reference": "/hpcnfs/scratch/DP/dfernand/Stamburri/PCGFS_PROJECT/ChIPseq/INPUT/{genome}/input_{control}.bam".format(
            genome = SAMPLES.GENOME[wildcards.sample], 
            control = wildcards.control) }
    return { "case": "02aln/{sample}.bam".format(sample=wildcards.sample),
             "reference": "/hpcnfs/scratch/DP/dfernand/Stamburri/PCGFS_PROJECT/ChIPseq/INPUT/{genome}/input_{control}.bam".format(
             genome = SAMPLES.GENOME[wildcards.sample], 
             control = wildcards.control),
             "spike": "02aln_dm/{sample}_spike.bam.clean".format(sample=wildcards.sample) }