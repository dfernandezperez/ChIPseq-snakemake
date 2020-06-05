def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()

rule cp_fastq_pe:
    input:
        get_fastq
    output:
        fastq1=temp("results/fastq/{sample}-{lane}.1.fastq.gz"),
        fastq2=temp("results/fastq/{sample}-{lane}.2.fastq.gz")
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
        temp("results/fastq/{sample}-{lane}.fastq.gz"),
    message:
        "Copying fastq files {input}"
    shell:
        """
        ln -s {input} {output}
        """


rule mergeFastq_pe:
    input:
        fw = lambda w: expand("results/fastq/{lane.sample}-{lane.lane}.1.fastq.gz", lane=units.loc[w.sample].itertuples()),
        rv = lambda w: expand("results/fastq/{lane.sample}-{lane.lane}.2.fastq.gz", lane=units.loc[w.sample].itertuples())
    output:
        fastq1 = temp("{tmp}/fastq/{{sample}}.1.fastq.gz".format(tmp=config["tmp"])),
        fastq2 = temp("{tmp}/fastq/{{sample}}.2.fastq.gz".format(tmp=config["tmp"]))
    log:
        "results/00log/fastp/{sample}.log"
    message:
        "Merging fastq files from {input}"
    shell:
        """
        cat {input.fw} > {output.fastq1}
        cat {input.rv} > {output.fastq2}
        """


rule mergeFastq_se:
    input:
        lambda w: expand("results/fastq/{lane.sample}-{lane.lane}.fastq.gz", lane=units.loc[w.sample].itertuples()),
    output:
        temp("{tmp}/fastq/{{sample}}.se.fastq.gz".format(tmp=config["tmp"]))
    log:
        "results/00log/fastp/{sample}.log"
    message:
        "Merging fastq files from {input}"
    shell:
        """
        cat {input} > {output}
        """


rule fastp_pe:
    input:
        fw = "{tmp}/fastq/{{sample}}.1.fastq.gz".format(tmp=config["tmp"]),
        rv = "{tmp}/fastq/{{sample}}.2.fastq.gz".format(tmp=config["tmp"])
    output:
        fastq1 = temp("{tmp}/fastq/trimmed/{{sample}}.1.fastq.gz".format(tmp=config["tmp"])),
        fastq2 = temp("{tmp}/fastq/trimmed/{{sample}}.2.fastq.gz".format(tmp=config["tmp"]))
    log:
        "00log/fastp/{sample}.log"
    threads:
        CLUSTER["fastp_pe"]["cpu"]
    params:
        fastp_params = config["params"]["fastp"]["pe"],
    message:
        "Processing fastq files from {input}"
    shadow:
        "minimal"
    shell:
        """
        fastp -i {input.fw} \
        -I {input.rv} \
        -o {output.fastq1} \
        -O {output.fastq2} \
        -w {threads} \
        {params.fastp_params} 2> {log}
        """


rule fastp_se:
    input:
        "{tmp}/fastq/{{sample}}.se.fastq.gz".format(tmp=config["tmp"])
    output:
        temp("{tmp}/fastq/trimmed/{{sample}}.se.fastq.gz".format(tmp=config["tmp"]))
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
    shell:
        """
        fastp -i {input} \
        -o {output} \
        -w {threads} {params.fastp_params} 2> {log}
        """