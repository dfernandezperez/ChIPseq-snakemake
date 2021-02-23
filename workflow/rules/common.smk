def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])


def is_spike(sample):
    return SAMPLES.SPIKE[sample] == True


# Get raw or trimmed reads based on trimming configuration
def get_fq(wildcards):
    if config["trimming"]:
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/trimmed{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)


def get_fq_spike(wildcards):
    if is_spike(**wildcards):
        if config["trimming"]:
            if not is_single_end(**wildcards):
                # paired-end sample
                return expand("{tmp}/fastq/trimmed{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
            # single end sample
            return "{tmp}/fastq/trimmed{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
        else:
            # no trimming, use raw reads
            if not is_single_end(**wildcards):
                # paired-end sample
                return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
            # single end sample
            return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)


# Get raw or trimmed reads based on trimming configuration. Used for fastqc
def get_fq_forward(wildcards):
    if config["trimming"]:
        if not is_single_end(wildcards.sample):
            # paired-end sample
            return "{tmp}/fastq/trimmed{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(wildcards.sample):
            # paired-end sample
            return "{tmp}/fastq/{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)


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
             "spike": "results/02aln_dm/{sample}_spike.bam.clean".format(sample=wildcards.sample),
             "ref_mm": "results/02aln/{input}.bam".format(input=wildcards.input),
             "ref_dm": "results/02aln_dm/{input}_spike.bam.clean".format(input=wildcards.input)
              }