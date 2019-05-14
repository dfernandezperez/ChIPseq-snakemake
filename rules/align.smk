rule align:
    input:
        get_trimmed
    output:
        bam = "02aln/{sample}.bam",
    threads:
        CLUSTER["align"]["cpu"]
    params:
        bowtie = "--chunkmbs 1024 -m 1 --best -S --no-unal -q " + config["idx_bt1_mm"],
        pe_params = "-I 10 -X 1000",
    message:
        "Aligning {input} with parameters {params.bowtie}"
    log:
       align   = "00log/alignments/{sample}.log",
       rm_dups = "00log/alignments/rm_dup/{sample}.log",
    benchmark:
        ".benchmarks/{sample}.align.benchmark.txt"
    run:
        n = len(input)
        if n == 1:
            reads = "{}".format(*input)
        else:
            reads = params.pe_params + " -1 {} -2 {}".format(*input)
        shell("""
        bowtie -p {threads} {params.bowtie} {reads} 2> {log.bowtie} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 5G -@ {threads} -T {output}.tmp -o {output} - 2>> {log.align}
        samtools index {output.} 2>> {log.align}
        """)


rule align_spike:
    input:
        get_trimmed_spike
    output:
        "02aln/{sample}_spike.bam"
    threads:
        CLUSTER["align"]["cpu"]
    params:
        bowtie = "--chunkmbs 1024 -m 1 --best -S --no-unal -q " + config["idx_bt1_spike"],
        pe_params = "-I 10 -X 1000",
    message:
        "Aligning {input} with parameters {params.bowtie}"
    log:
       align   = "00log/alignments/{sample}_spike.log",
       rm_dups = "00log/alignments/rm_dup/{sample}_spike.log",
    benchmark:
        ".benchmarks/{sample}.alignSpike.benchmark.txt"
    run:
        n = len(input)
        if n == 1:
            reads = "{}".format(*input)
        else:
            reads = params.pe_params + " -1 {} -2 {}".format(*input)
        shell("""
        bowtie -p {threads} {params.bowtie} {reads} 2> {log.bowtie} \
        | samblaster --removeDups 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 5G -@ {threads} -T {output}.tmp -o {output} - 2>> {log.align}
        samtools index {output} 2>> {log.align}
        """)


rule clean_spike:
    input:
        mm = "02aln/{sample}.bam",
        spike = "02aln_dm/{sample}_spike.bam"
    output:
        mm = "02aln/{sample}.clean.bam",
        dm = "02aln_dm/{sample}_spike.clean.bam"
    log:
        "00log/alignments/{sample}.removeSpikeDups"
    shell:
        """
        python scripts/remove_spikeDups.py {input}       
        mv {input.mm}.clean {output.mm}; mv {input.dm}.clean {output.dm}
        samtools index {output.mm}; samtools index {output.dm}
        """


rule bam2bigwig:
    input: 
        unpack(get_bam_cntrl)
    output:  
        "06bigwig/{sample}_{control}-input.bw"
    params: 
        read_exten = config['read_extension']
    log: 
        "00log/{sample}_{control}-input_bigwig.bam2bw"
    threads: 
        CLUSTER["bam2bigwig"]["cpu"]
    message: 
        "making input subtracted bigwig for sample {wildcards.sample} with input {input.reference}"
    benchmark:
        ".benchmarks/{sample}_{control}.bam2bw.benchmark.txt"
    shell:
        """
        python scripts/bam2bigwig.py --case {input.case} \
        --reference {input.reference} \
        --bigwig {output} \
        --extReads {params.read_exten} \
        --threads {threads} 2> {log}
        """

rule bam2bigwig_spike:
    input: 
        unpack(get_bam_cntrl)
    output:
        "06bigwig/{sample}_{control}-inputSpike.bw"
    params:
        chr_sizes  = config["genome"]["chr_sizes"],
        read_exten = config['read_extension']
    log: 
        "00log/{sample}_{control}-input_bigwig.bam2bw"
    threads:
        CLUSTER["bam2bigwig"]["cpu"]
    message:
        "making spike-normalized input subtracted bigwig for sample {wildcards.sample} with input {input.reference}"
    benchmark:
        ".benchmarks/{sample}_{control}.bam2bw-spike.benchmark.txt"
    shell:
        """
        python scripts/bam2bigwig_spike.py --case {input.case} \
        --reference {input.reference} \
        --spike {input.dm} \
        --bigwig {output} \
        --extReads {params.read_exten} \
        --chrSizes {params.chr_sizes} \
        --threads {threads} 2> {log}
        """

# rule bigwig2server:
#     input: 
#         bw       = "06bigwig/{sample}_{control}-input.bw",
#         flagstat = "02aln/{sample}.bam.flagstat"
#     output:
#         temp("temp_file_{sample}_{control}.txt")
#     params:
#         user     = lambda wildcards : SAMPLES.USER[wildcards.sample],
#         antibody = lambda wildcards : SAMPLES.AB[wildcards.sample],
#         genome   = lambda wildcards : SAMPLES.GENOME[wildcards.sample],
#         run      = lambda wildcards : SAMPLES.RUN[wildcards.sample],
#         chip     = lambda wildcards : str("ChIPseq") if SAMPLES.SPIKE[wildcards.sample] == False else str("ChIPseqSpike")
#     run:
#         with open (input.flagstat, "r") as f:
#             line         = f.readlines()[4] # fifth line contains the number of mapped reads
#             match_number = re.match(r'(\d.+) \+.+', line)
#             total_reads  = int(match_number.group(1))
#         shell("cp {input} \
#             /hpcnfs/data/DP/UCSC_tracks/Data/bigWig/{sample}_{control}_{user}_{nreads}_{chip}_{antibody}_{genome}_{run}.bigWig".format(
#             input    = input.bw,
#             sample   = wildcards.sample,
#             control  = wildcards.control,
#             user     = params.user,
#             nreads   = total_reads,
#             chip     = params.chip,
#             antibody = params.antibody,
#             genome   = params.genome,
#             run      = params.run))
#         shell("touch {output}".format(output = output))