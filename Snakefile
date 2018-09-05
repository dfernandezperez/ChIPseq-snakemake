import pandas as pd
singularity: "/hpcnfs/data/DP/Singularity/dfernandezperez-ChIPseq-software-master-latest.simg"

#######################################################################################################################
### Load sample sheet and cluster configuration, config file
#######################################################################################################################
configfile: "config.yaml"
CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = pd.read_table(config['samples']).set_index("NAME", drop=False)

#######################################################################################################################
### Process a little bit the input files. Check if there are samples with spike-in and define all the target output files
#######################################################################################################################
ALL_SAMPLES = FILES.NAME
# create 2 lists, one containing the samples with spike in and the other without spike-in
NOSPIKE_SAMPLES = list( FILES[FILES.SPIKE == False].NAME )
SPIKE_SAMPLES = list( FILES[FILES.SPIKE == True].NAME )

# Now, create the path of the fastq files that will be used as input by bowtie for all samples
# Finally, create a dictionary that has as keys the sampleName and as values the fastq paths.
#This way then we can acess to each fastq by using the sample name
SPIKE_FASTQ = expand("fastq/{sample}.fastq", sample = SPIKE_SAMPLES)
NOSPIKE_FASTQ = expand("fastq/{sample}.fastq", sample = NOSPIKE_SAMPLES)
SPIKE_FASTQ = dict(zip(SPIKE_SAMPLES, SPIKE_FASTQ))
NOSPIKE_FASTQ = dict(zip(NOSPIKE_SAMPLES, NOSPIKE_FASTQ))

# Do the same but for the bam files, that will be used for the bigwig creation (spike in have to be normalized in a diff way)
SPIKE_BAM = expand("02aln/{sample}.bam", sample = SPIKE_SAMPLES)
NOSPIKE_BAM = expand("02aln/{sample}.bam", sample = NOSPIKE_SAMPLES)
SPIKE_BAM = dict(zip(SPIKE_SAMPLES, SPIKE_BAM))
NOSPIKE_BAM = dict(zip(NOSPIKE_SAMPLES, NOSPIKE_BAM))


#######################################################################################################################
### DETERMINE ALL THE OUTPUT FILES TO RUN SNAKEMAKE
#######################################################################################################################
# Get the names of the input files for calling peak and creating bigwig files
CONTROLS = FILES.INPUT

ALL_FASTQC = expand("01fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES)
ALL_BAM = expand("02aln/{sample}.bam", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("02aln/{sample}.bam.flagstat", sample = ALL_SAMPLES)
ALL_PEAKS = expand("03peak_macs2/{sample}_{control}-input/{sample}_peaks.narrowPeak", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_PHANTOM = expand("04phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
ALL_PEAKANNOT = expand("07peak_annot/{sample}_{control}-input/{sample}_peaks_p10.annot", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_BIGWIG = expand( "05bigwig/{sample}_{control}-input.bw", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_GCBIAS = expand("06gcBias/{sample}_{control}-input_GCbias.pdf", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_QC = ["10multiQC/multiQC_log.html"]

#######################################################################################################################
### DEFINE LOCAL RULES TO RUN THE WHOLE PIPELINE OR JUST A SUBSET OF IT
#######################################################################################################################
localrules: all
rule all:
    input: ALL_FASTQC + ALL_BAM + ALL_FLAGSTAT + ALL_PEAKS + ALL_PHANTOM + ALL_BIGWIG + ALL_PEAKANNOT + ALL_QC + ALL_GCBIAS
    
#######################################################################################################################
### FASTQC, ALIGNMENT + DEDUP + SABM2BAM SORTED, INDEX BAM
#######################################################################################################################
## get a list of fastq.gz files for each sample
def get_fastq(wildcards):
    return FILES.FASTQ[wildcards.sample]

rule merge_fastqs:
    input: 
        get_fastq
    output: 
        temp("fastq/{sample}.fastq")
    log: 
        "00log/{sample}_unzip"
    threads: 
        CLUSTER["merge_fastqs"]["cpu"]
    params: 
        fastp = "/hpcnfs/data/DP/software/fastp", 
        fastp_params = "-j 00log/{sample}_fastp.json -h 00log/{sample}_fastp.html --stdin -t 1 -A -Q -L "
    message: 
        "merging fastq files {input}/*.fastq.gz into {output}"
    shell:
        """
        zcat {input}/*fastq.gz | {params.fastp} -o {output} -w {threads} {params.fastp_params} 2> {log}
        rm 00log/{wildcards.sample}_fastp.json 00log/{wildcards.sample}_fastp.html
        """

rule fastqc:
    input:  
        "fastq/{sample}.fastq"
    output: 
        "01fqc/{sample}_fastqc.zip"
    log:    
        "00log/{sample}_fastqc"
    threads: 
        CLUSTER["fastqc"]["cpu"]
    message: 
        "Running fastqc for {input}"
    shell:
        """
        fastqc -o 01fqc -f fastq -t {threads} --noextract {input} 2> {log}
        """

# Check if there's spike in, get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and duplicated reads by samblaster -r
# samblaster should run before samtools sort
if NOSPIKE_SAMPLES:
    rule align:
        input:
            lambda wildcards: NOSPIKE_FASTQ[wildcards.sample]
        output:
            temp("02aln/{sample}.bam")
        threads:
            CLUSTER["align"]["cpu"]
        params:
            bowtie = "--chunkmbs 1024 -m 1 --best -S " + config["idx_bt1_mm"]
        message:
            "Aligning {input} with parameters {params.bowtie}"
        log:
            bowtie = "00log/{sample}.align",
            markdup = "00log/{sample}.markdup",
            sort_index = "00log/{sample}.sortIndex_bam"
        shell:
            """
            bowtie -p {threads} {params.bowtie} -q {input} 2> {log.bowtie} \
            | samblaster --removeDups 2> {log.markdup} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m 2G -@ {threads} -T {output}.tmp -o {output} - 2> {log.sort_index}
            samtools index {output} 2>> {log.sort_index}
            """

if SPIKE_SAMPLES:
    rule align_spike:
        input:
            lambda wildcards: SPIKE_FASTQ[wildcards.sample]
        output:
            mm = temp("02aln/{sample}.bam"),
            dm = "02aln_dm/{sample}_dm.bam"
        threads:
            CLUSTER["align"]["cpu"]
        params:
            bowtie_mm = "--chunkmbs 1024 -m 1 --best -S " + config["idx_bt1_mm"],
            bowtie_dm = "--chunkmbs 1024 -m 1 --best -S " + config["idx_bt1_dm"]
        message:
            "Aligning {input} with parameters {params.bowtie_mm}"
        log:
            bowtie = ["00log/{sample}.align", "00log/{sample}.dm_align"],
            markdup = "00log/{sample}.markdup",
            sort_index = "00log/{sample}.sortIndex_bam"
        shell:
            """
            bowtie -p {threads} {params.bowtie_mm} -q {input} 2> {log.bowtie[0]} \
            | samblaster --removeDups 2> {log.markdup} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m 2G -@ {threads} -T {output.mm}.tmp -o {output.mm} - 2> {log.sort_index}
            samtools index {output.mm} 2>> {log.sort_index}

            bowtie -p {threads} {params.bowtie_dm} -q {input} 2> {log.bowtie[1]} \
            | samblaster --removeDups 2> {log.markdup} \
            | samtools view -Sb -F 4 - \
            | samtools sort -m 2G -@ {threads} -T {output.dm}.tmp -o {output.dm} - 2> {log.sort_index}
            samtools index {output.dm} 2>> {log.sort_index}

            python scripts/remove_spikeDups.py {output.mm} {output.dm}
            
            mv {output.mm}.clean {output.mm}; mv {output.dm}.clean {output.dm}

            samtools index {output.mm}; samtools index {output.dm}
            """

# This one is to have the number of mapped reads after removing the duplicates in a format that multiQC can recognize
rule flagstat_bam:
    input:  
        "02aln/{sample}.bam"
    output: 
        "02aln/{sample}.bam.flagstat"
    log:    
        "00log/{sample}.flagstat_bam"
    message: 
        "flagstat_bam {input}"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

#######################################################################################################################
### PEAK CALLING WITH MACS2 AND PEAK QUALITY STATS WITH SPP
#######################################################################################################################
rule call_peaks:
    input: 
        case = "02aln/{sample}.bam", 
        reference = "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/input_{control}.bam"
    output: 
        narrowPeak = "03peak_macs2/{sample}_{control}-input/{sample}_peaks.narrowPeak",
        bed_p10 = "03peak_macs2/{sample}_{control}-input/{sample}_peaks_p10.bed",
        bdg = temp("03peak_macs2/{sample}_{control}-input/{sample}_treat_pileup.bdg")
    log:
        "00log/{sample}_{control}-input_macs2.log"
    params:
        out_dir = "03peak_macs2/{sample}_{control}-input",
        macs2_params = "--keep-dup all --format BAM --nomodel --bdg"
    message: 
        "call_peaks macs2 with input {input.reference} for sample {input.case}"
    shell:
        """
        macs2 callpeak --treatment {input.case} \
            --control {input.reference} \
            --gsize {config[macs2_gsize]} \
            --outdir {params.out_dir} \
            --name {wildcards.sample} \
            --pvalue {config[macs2_pvalue]} \
            {params.macs2_params} &> {log}            
        awk "\$8>=10" {output.narrowPeak} | cut -f1-4,8 > {output.bed_p10}
        rm *_control_lambda.bdg # Remove useless bedgraph created by macs2
        """

rule phantom_peak_qual:
    input: 
        "02aln/{sample}.bam"
    output:
        "04phantompeakqual/{sample}_phantom.txt"
    log:
        "00log/{sample}_phantompeakqual.log"
    threads:
        CLUSTER["phantom_peak_qual"]["cpu"]
    params:
        out_dir = "04phantompeakqual"
    message:
        "Running phantompeakqual for {wildcards.sample}"
    shell:
        """
        /usr/local/bin/Rscript  scripts/run_spp_nodups.R \
        -c={input[0]} -savp -rf -p={threads} -odir={params.out_dir}  -out={output} -tmpdir={params.out_dir}  2> {log}
        """

rule peakAnnot:
    input :
        "03peak_macs2/{sample}_{control}-input/{sample}_peaks_p10.bed"
    output:
        annot = "07peak_annot/{sample}_{control}-input/{sample}_peaks_p10.annot",
        promo_bed_targets = "07peak_annot/{sample}_{control}-input/{sample}_peaks_p10_promoTargets.bed",
        promoTargets = "07peak_annot/{sample}_{control}-input/{sample}_peaks_p10_promoTargets.txt",
        promoBed = "07peak_annot/{sample}_{control}-input/{sample}_peaks_p10_promoPeaks.bed",
        distalBed = "07peak_annot/{sample}_{control}-input/{sample}_peaks_p10_distalPeaks.bed"
    log: 
        "00log/{sample}_{control}-input_peakanot"
    message:
        "Annotating peaks for {wildcards.sample}"
    shell:
        """
        /usr/local/bin/Rscript scripts/PeakAnnot.R {input} {config[promoter][bTSS]} {config[promoter][aTSS]}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} {output.distalBed}
        """

#######################################################################################################################
### BAM TO BIGWIG WITH DEEPTOOLS, GC BIAS
#######################################################################################################################
if NOSPIKE_SAMPLES:
    rule bam2bigwig:
        input: 
            case = lambda wildcards: NOSPIKE_BAM[wildcards.sample],
            reference = "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/input_{control}.bam"
        output:  
            "05bigwig/{sample}_{control}-input.bw"
        params: 
            read_exten = config['read_extension']
        log: 
            "00log/{sample}_{control}-input_bigwig.bam2bw"
        threads: 
            CLUSTER["bam2bigwig"]["cpu"]
        message: 
            "making input subtracted bigwig for sample {wildcards.sample} with input {input.reference}"
        shell:
            """
            python scripts/bam2bigwig.py --case {input.case} \
            --reference {input.reference} \
            --bigwig {output} \
            --extReads {params.read_exten} \
            --threads {threads} &> {log}
            """

if SPIKE_SAMPLES:
    rule bam2bigwig_spike:
        input: 
            case = lambda wildcards: SPIKE_BAM[wildcards.sample],
            dm = "02aln_dm/{sample}_dm.bam",
            reference = "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/input_{control}.bam"
        output:
            "05bigwig/{sample}_{control}-input.bw"
        params:
            chr_sizes = config["chr_sizes"],
            read_exten = config['read_extension']
        log: 
            "00log/{sample}_{control}-input_bigwig.bam2bw"
        threads:
            CLUSTER["bam2bigwig"]["cpu"]
        message:
            "making spike-normalized input subtracted bigwig for sample {wildcards.sample} with input {input.reference}"
        shell:
            """
            python scripts/bam2bigwig_spike.py --case {input.case} \
            --reference {input.reference} \
            --spike {input.dm} \
            --bigwig {output} \
            --extReads {params.read_exten} \
            --chrSizes {params.chr_sizes} \
            --threads {threads} &> {log}
            """

rule GC_bias:
    input: 
        bam = "02aln/{sample}.bam",
        bed = "03peak_macs2/{sample}_{control}-input/{sample}_peaks_p10.bed"
    output: 
        pdf = "06gcBias/{sample}_{control}-input_GCbias.pdf",
        freq_txt = "06gcBias/{sample}_{control}-input_GCbias.txt"
    log:
        "00log/{sample}_{control}-input_GCbias.log"
    params:
        repeatMasker = config['rep_masker'],
        tempBed = "06gcBias/{sample}_Repeatmasker.bed.tmp",
        bit_file = config["2bit"],
        egenome_size = config["egenome_size"]
    threads:
        CLUSTER["GC_bias"]["cpu"]
    message:
        "Computing GC bias for sample {wildcards.sample}"
    shell:
        """
        bedops -u {input.bed} {params.repeatMasker} > {params.tempBed}
        computeGCBias -b {input.bam} \
            -p {threads} \
            --effectiveGenomeSize {params.egenome_size} \
            -g {params.bit_file} \
            -l 200 \
            -bl {params.tempBed} \
            --biasPlot {output.pdf} \
            --GCbiasFrequenciesFile {output.freq_txt} 2> {log}
        rm -f {params.tempBed}
        """

#######################################################################################################################
### FINAL HTML REPORT
#######################################################################################################################
rule multiQC:
    input:
        expand("00log/{sample}.align", sample = ALL_SAMPLES),
        expand("01fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES),
        expand("02aln/{sample}.bam.flagstat", sample = ALL_SAMPLES)
    output: 
        "10multiQC/multiQC_log.html"
    params:
        log_name = "multiQC_log"
    log:
        "00log/multiqc.log"
    message:
        "multiqc for all logs"
    shell:
        """
        multiqc {input} -o 10multiQC -f -v -n {params.log_name} 2> {log}
        """