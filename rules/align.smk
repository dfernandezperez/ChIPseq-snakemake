def set_reads(wildcards, input):
        n = len(input)
        if n == 1:
            reads = "{}".format(*input)
            return reads
        else:
            reads = config["params"]["bowtie"]["pe"] + " -1 {} -2 {}".format(*input)
            return reads


rule align:
    input:
        get_trimmed
    output:
         bam   = temp("02aln/{sample}.bam.tmp"),
         index = temp("02aln/{sample}.bam.tmp.bai")
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index  = config["ref"]["index"],
        bowtie = config["params"]["bowtie"]["global"],
        reads  = set_reads,
    message:
        "Aligning {input} with parameters {params.bowtie}"
    log:
       align   = "00log/alignments/{sample}.log",
       rm_dups = "00log/alignments/rm_dup/{sample}.log",
    benchmark:
        ".benchmarks/{sample}.align.benchmark.txt"
    shell:
        """
        bowtie -p {threads} {params.bowtie} {params.index} {params.reads} 2> {log.align} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 5G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam}
        """


rule align_spike:
    input:
        get_trimmed_spike
    output:
        bam   = temp("02aln_dm/{sample}_spike.bam"),
        index = temp("02aln_dm/{sample}_spike.bam.bai")
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index  = config["ref"]["index_spike"],
        bowtie = config["params"]["bowtie"]["global"],
        reads  = set_reads,
    message:
        "Aligning {input} with parameters {params.bowtie}"
    log:
       align   = "00log/alignments/{sample}_spike.log",
       rm_dups = "00log/alignments/rm_dup/{sample}_spike.log",
    benchmark:
        ".benchmarks/{sample}.alignSpike.benchmark.txt"
    shell:
        """
        bowtie -p {threads} {params.bowtie} {params.index} {params.reads} 2> {log.align} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 5G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam}
        """


rule clean_spike:
    input:
        mm          = "02aln/{sample}.bam.tmp",
        spike       = "02aln_dm/{sample}_spike.bam",
        mm_index    = "02aln/{sample}.bam.tmp.bai",
        spike_index = "02aln_dm/{sample}_spike.bam.bai",
    output:
        mm    = temp("02aln/{sample}.bam.tmp.clean"),
        spike = "02aln_dm/{sample}_spike.bam.clean"
    log:
        "00log/alignments/{sample}.removeSpikeDups"
    shell:
        """
        python scripts/remove_spikeDups.py {input} &> {log}      
        mv {input.mm}.temporary {output.mm}; mv {input.spike}.temporary {output.spike}
        samtools index {output.spike}
        """

# Dummy rule to change the name of the bam files to be able to 
# have the same name structure in spike-in and non-spiked samples
rule update_bam:
    input:
        get_bam
    output:
        "02aln/{sample}.bam",
    log:
        "00log/alignments/{sample}.update_bam"
    shell:
        """
        cp {input} {output}
        samtools index {output} 2>> {log}
        """


def set_reads_spike(wildcards, input):
        n = len(input)
        assert n == 2 or n == 3, "input->sample must have 2 (sample + input) or 3 (sample + input + spike) elements"
        if n == 2:
            reads = "scripts/bam2bigwig.py"
            return reads
        if n == 3:
            reads = "scripts/bam2bigwig_spike.py --spike {} --chrSizes ".format(input.spike) + config["ref"]["chr_sizes"]
            return reads


rule bam2bigwig:
    input: 
        unpack(get_bam_cntrl)
    output:  
        "06bigwig/{sample}_{control}-input.bw"
    params: 
        read_exten = set_read_extension,
        reads = set_reads_spike,
    log: 
        "00log/bam2bw/{sample}_{control}-input_bigwig.bam2bw"
    threads: 
        CLUSTER["bam2bigwig"]["cpu"]
    message: 
        "making input subtracted bigwig for sample {wildcards.sample} with input {input.reference}"
    benchmark:
        ".benchmarks/{sample}_{control}.bam2bw.benchmark.txt"
    shell:
        """
        python {params.reads} \
        --case {input.case} \
        --reference {input.reference} \
        --bigwig {output} \
        --threads {threads} &> {log}
        """


def set_reads_spike2(wildcards, input):
        n = len(input)
        assert n == 1 or n == 2, "input->sample must have 1 (sample) or 2 (sample + spike) elements"
        if n == 1:
            reads = "scripts/bam2bigwig_noSubtract.py"
            return reads
        if n == 2:
            reads = "scripts/bam2bigwig_spike_noSubtract.py --spike {} --chrSizes ".format(input.spike) + config["ref"]["chr_sizes"]
            return reads

def input_noSubstract(w):
    if not is_spike(w.sample):
        return { "case": "02aln/{sample}.bam".format(sample=w.sample) }
    return { "case": "02aln/{sample}.bam".format(sample=w.sample),
             "spike": "02aln_dm/{sample}_spike.bam.clean".format(sample=w.sample) }

rule bam2bigwig_noSubstract:
    input: 
        unpack(input_noSubstract)
    output:  
        "06bigwig/noSubtract/{sample}.bw"
    params: 
        read_exten = set_read_extension,
        reads = set_reads_spike2,
    log: 
        "00log/bam2bw/{sample}_bigwig.bam2bw"
    threads: 
        CLUSTER["bam2bigwig"]["cpu"]
    message: 
        "making input subtracted bigwig for sample {wildcards.sample}"
    shell:
        """
        python {params.reads} \
        --case {input.case} \
        --bigwig {output} \
        --threads {threads} &> {log}
        """

rule bigwig2server:
    input: 
        bw       =  "06bigwig/noSubtract/{sample}.bw",
        samblaster = "00log/alignments/rm_dup/{sample}.log"
    output:
        temp("temp_file_{sample}_{control}.txt")
    params:
        user     = lambda wildcards : SAMPLES.USER[wildcards.sample],
        antibody = lambda wildcards : SAMPLES.AB[wildcards.sample],
        genome   = lambda wildcards : SAMPLES.GENOME[wildcards.sample],
        run      = lambda wildcards : SAMPLES.RUN[wildcards.sample],
        chip     = lambda wildcards : str("ChIPseq") if SAMPLES.SPIKE[wildcards.sample] == False else str("ChIPseqSpike")
    run:
        with open (input.samblaster, "r") as f:
            line         = f.readlines()[4] # fifth line contains the number of mapped reads
            match_number = re.match(r'samblaster: Removed (\d.+) of (\d.+) \(', line)
            total_reads  = int(match_number.group(2))-int(match_number.group(1))
        shell("cp {input} \
            /hpcnfs/data/DP/UCSC_tracks/Data/bigWig/{sample}_{control}_{user}_{nreads}_{chip}_{antibody}_{genome}_{run}.bigWig".format(
            input    = input.bw,
            sample   = wildcards.sample,
            control  = wildcards.control,
            user     = params.user,
            nreads   = total_reads,
            chip     = params.chip,
            antibody = params.antibody,
            genome   = params.genome,
            run      = params.run))
        shell("touch {output}".format(output = output))