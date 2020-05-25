def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])


def is_spike(sample):
    return SAMPLES.SPIKE[sample] == True


def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], tmp=config["tmp"], **wildcards)
    # single end sample
    return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp=config["tmp"], **wildcards)


def get_trimmed_spike(wildcards):
    if is_spike(**wildcards):
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], tmp=config["tmp"], **wildcards)
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp=config["tmp"], **wildcards)


# Used to use just the forward reads to create the fastqc report
def get_trimmed_forward(wildcards):
    if not is_single_end(**wildcards):
        # paired-end sample
        return "{tmp}/fastq/{sample}.1.fastq.gz".format(tmp=config["tmp"], **wildcards)
    # single end sample
    return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp=config["tmp"], **wildcards)


def get_bam(wildcards):
    if not is_spike(**wildcards):
        return "results/02aln/{sample}.bam.tmp".format(**wildcards)
    return "results/02aln/{sample}.bam.tmp.clean".format(**wildcards)


def set_read_extension(wildcards):
    if is_single_end(wildcards.sample):
        return "--extendReads " + str(config['bam2bigwig']['read_extension'])
    return "--extendReads"


def get_bam_cntrl(wildcards):
    if not is_spike(wildcards.sample):
        return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample),
            "reference": "results/02aln/{control}.bam" }
    return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample),
             "reference": "results/02aln/{control}.bam",
             "spike": "results/02aln_dm/{sample}_spike.bam.clean".format(sample=wildcards.sample) }


def get_bam_spike(wildcards):
    if not is_spike(wildcards.sample):
        return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample) }
    return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample),
             "spike": "results/02aln_dm/{sample}_spike.bam.clean".format(sample=wildcards.sample) }