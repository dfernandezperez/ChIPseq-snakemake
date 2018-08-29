#######################################################################################################################
### Load sample sheet and cluster configuration, config file
#######################################################################################################################
configfile: "config.yaml"
CLUSTER = json.load(open(config['CLUSTER_JSON']))
FILES = json.load(open(config['SAMPLES_JSON']))
ALL_SAMPLES = sorted(FILES.keys())


#######################################################################################################################
### Process a little bit the input files. Check if there are samples with spike-in and define all the target output files
#######################################################################################################################
# create 2 lists, one containing the samples with spike in and the other without spike-in
NOSPIKE_SAMPLES = []
SPIKE_SAMPLES = []
for sample in ALL_SAMPLES:
    if "TRUE" == FILES[sample]["SPIKE"]: # If the sample has spike in, the field FILES[sample]["SPIKE"] must contain "TRUE"
        SPIKE_SAMPLES.append(sample)
    elif "FALSE" == FILES[sample]["SPIKE"]:
        NOSPIKE_SAMPLES.append(sample)

# Now, create the path of the fastq files that will be used as input by bowtie for all samples
SPIKE_FASTQ = expand("fastq/{sample}.fastq", sample = SPIKE_SAMPLES)
NOSPIKE_FASTQ = expand("fastq/{sample}.fastq", sample = NOSPIKE_SAMPLES)

# Finally, create a dictionary that has as keys the sampleName and as values the fastq paths.
#This way then we can acess to each fastq by using the sample name
SPIKE_FASTQ = dict(zip(SPIKE_SAMPLES, SPIKE_FASTQ))
NOSPIKE_FASTQ = dict(zip(NOSPIKE_SAMPLES, NOSPIKE_FASTQ))

# Do the same but for the bam files, that will be used for the bigwig creation (spike in have to be normalized in a diff way)
SPIKE_BAM = expand("02aln/{sample}.sorted.bam", sample = SPIKE_SAMPLES)
NOSPIKE_BAM = expand("02aln/{sample}.sorted.bam", sample = NOSPIKE_SAMPLES)
SPIKE_BAM = dict(zip(SPIKE_SAMPLES, SPIKE_BAM))
NOSPIKE_BAM = dict(zip(NOSPIKE_SAMPLES, NOSPIKE_BAM))



#######################################################################################################################
### DETERMINE ALL THE OUTPUT FILES TO RUN SNAKEMAKE
#######################################################################################################################
# Get the names of the input files for calling peak and creating bigwig files
CONTROLS = [FILES[x]["INPUT"] for x in ALL_SAMPLES]


ALL_FASTQC = expand("01fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES)
ALL_BAM = expand("02aln/{sample}.sorted.bam", sample = ALL_SAMPLES)
ALL_FLAGSTAT = expand("02aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
ALL_PEAKS = expand("03peak_macs2/{sample}_{control}-input/{sample}_peaks.xls", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_PHANTOM = expand("04phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
ALL_PEAKANNOT = expand("peak_annot/{sample}_{control}-input/{sample}_peaks_p10.annot", zip, sample=ALL_SAMPLES, control=CONTROLS)
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
    return FILES[wildcards.sample]["FASTQ"]

rule merge_fastqs:
    input: get_fastq
    output: temp("fastq/{sample}.fastq")
    log: "00log/{sample}_unzip"
    threads: CLUSTER["merge_fastqs"]["cpu"]
    params: fastp = " /hpcnfs/data/DP/software/fastp"
    message: "merging fastq files {input}/*.fastq.gz into {output}"
    shell:
        """
        set +u; source activate DPbase; set -u
        zcat {input}/*fastq.gz | {params.fastp} --stdin -o {output} -w {threads} -t 1 -A -Q -L 2> {log}
        rm fastp.html fastp.json # Remove the stupid report files created by fastp
        """

rule fastqc:
    input:  "fastq/{sample}.fastq"
    output: "01fqc/{sample}_fastqc.zip"
    log:    "00log/{sample}_fastqc"
    threads: CLUSTER["fastqc"]["cpu"]
    message: "Running fastqc for {input}"
    shell:
        """
        set +u; source activate DPbase; set -u
        fastqc -o 01fqc -f fastq -t {threads} --noextract {input} 2> {log}
        """

# Check if there's spike in, get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and duplicated reads by samblaster -r
# samblaster should run before samtools sort
if  SPIKE_SAMPLES:
    rule align_spike:
        input: lambda wildcards: SPIKE_FASTQ[wildcards.sample]
        output: mm = temp("02aln/{sample}.bam"), dm = "02aln_dm6/{sample}_dm6.bam", "00log/{sample}.align"
        threads: CLUSTER["align"]["cpu"]
        params:
                bowtie_mm = "--chunkmbs 1024 -m 1 --best " + config["idx_bt1_mm"],
                bowtie_dm = "--chunkmbs 1024 -m 1 --best " + config["idx_bt1_dm"]
        message: "Aligning {input} with parameters {params.bowtie_mm}"
        log:
            bowtie = ["00log/{sample}.align", "00log/{sample}.dm6_align"],
            markdup = "00log/{sample}.markdup"
        shell:
            """
            set +u; source activate DPbase; set -u
            bowtie -p {threads} {params.bowtie_mm} -q {input} -S 2> {log.bowtie[0]} \
            | samblaster --removeDups 2> {log.markdup} \
            | samtools view -Sb -F 4 - > {output.mm}.tmp

            bowtie -p {threads} {params.bowtie_dm} -q {input} -S 2> {log.bowtie[1]} \
            | samblaster --removeDups 2> {log.markdup} \
            | samtools view -Sb -F 4 - > {output.dm}.tmp

            python scripts/remove_spikeDups.py {output[0]} {output[1]}
            
            mv {output.mm}.clean {output.mm}; mv {output.dm}.clean {output.dm}
            """

if NOSPIKE_SAMPLES:
    rule align:
        input: lambda wildcards: NOSPIKE_FASTQ[wildcards.sample]
        output: "02aln/{sample}.bam", "00log/{sample}.align"
        threads: CLUSTER["align"]["cpu"]
        params:
                bowtie = "--chunkmbs 1024 -m 1 --best " + config["idx_bt1_mm"]
        message: "Aligning {input} with parameters {params.bowtie}"
        log:
            bowtie = "00log/{sample}.align",
            markdup = "00log/{sample}.markdup"
        shell:
            """
            set +u; source activate DPbase; set -u
            bowtie -p {threads} {params.bowtie} -q {input} -S 2> {log.bowtie} \
            | samblaster --removeDups 2> {log.markdup} \
            | samtools view -Sb -F 4 - > {output[0]}
            """

rule sortIndex_bam:
    input:  "02aln/{sample}.bam"
    output: "02aln/{sample}.sorted.bam"
    log:    "00log/{sample}.sortIndex_bam"
    threads: CLUSTER["sortIndex_bam"]["cpu"]
    message: "sorting and indexing {input}"
    shell:
        """
        set +u; source activate DPbase; set -u
        samtools sort -m 2G -@ {threads} -T {output}.tmp -o {output} {input} 2> {log}
        samtools index {output} 2>> {log}
        """

# This one is to have the number of mapped reads after removing the duplicates in a format that multiQC can recognize
rule flagstat_bam:
    input:  "02aln/{sample}.sorted.bam"
    output: "02aln/{sample}.sorted.bam.flagstat"
    log:    "00log/{sample}.flagstat_bam"
    message: "flagstat_bam {input}"
    shell:
        """
        set +u; source activate DPbase; set -u
        samtools flagstat {input} > {output} 2> {log}
        """

#######################################################################################################################
### PEAK CALLING WITH MACS2 AND PEAK QUALITY STATS WITH SPP
#######################################################################################################################
rule call_peaks:
    input: case="02aln/{sample}.sorted.bam", reference = "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/input_{control}.bam"
    output: peaks_xls="03peak_macs2/{sample}_{control}-input/{sample}_peaks.xls",
            narrowPeak="03peak_macs2/{sample}_{control}-input/{sample}_peaks.narrowPeak",
            bed_p10="03peak_macs2/{sample}_{control}-input/{sample}_peaks_p10.bed"
    log: "00log/{sample}_{control}-input_macs2.log"
    params:
            jobname = "{sample}", prefix = "03peak_macs2/{sample}_{control}-input"
    message: "call_peaks macs2 with input {input.reference} for sample {input.case}"
    shell:
        """
        set +u; source activate DPbase; set -u
        macs2 callpeak -t {input.case} \
            -c {input.reference} --keep-dup all -f BAM -g {config[macs2_g]} \
            --outdir {params.prefix} -n {wildcards.sample} -p {config[macs2_pvalue]} --nomodel &> {log}
        awk '$8>=10' {output.narrowPeak} | cut -f1-4,8 > {output.bed_p10}
        """

rule phantom_peak_qual:
    input: "02aln/{sample}.sorted.bam"
    output: "04phantompeakqual/{sample}_phantom.txt"
    log: "00log/{sample}_phantompeakqual.log"
    threads: CLUSTER["phantom_peak_qual"]["cpu"]
    params: jobname = "{sample}", prefix = "04phantompeakqual"
    message: "phantompeakqual for {params.jobname}"
    shell:
        """
        set +u; source activate DPbase; set -u
        Rscript  /hpcnfs/scratch/DP/dfernand/phantompeakqualtools/run_spp_nodups.R \
        -c={input[0]} -savp -rf -p={threads} -odir={params.prefix}  -out={output} -tmpdir={params.prefix}  2> {log}
        """

rule peakAnnot:
    input : rules.call_peaks.output.bed_p10
    output: annot="peak_annot/{sample}_{control}-input/{sample}_peaks_p10.annot",
            promo_bed_targets="peak_annot/{sample}_{control}-input/{sample}_peaks_p10_promoTargets.bed",
            promoTargets="peak_annot/{sample}_{control}-input/{sample}_peaks_p10_promoTargets.txt",
            promoBed="peak_annot/{sample}_{control}-input/{sample}_peaks_p10_promoPeaks.bed",
            distalBed="peak_annot/{sample}_{control}-input/{sample}_peaks_p10_distalPeaks.bed"
    log: "00log/{sample}_{control}-input_peakanot"
    message: "Annotating peaks for {wildcards.sample}"
    shell:
        """
        set +u; source activate DPbase; set -u
        /hpcnfs/scratch/DP/dfernand/Scripts/PeakAnnot.R {input} {output.annot} 2500 2500 2> {log}
        awk '$8 == "Promoter" {{OFS="\t"; if($13 == 1) $13="+"; else $13="-" ; print $1,$10,$11,$15,$18,$13}}' {output.annot} | sort -u > {output.promo_bed_targets}
        awk '$8 == "Promoter"' {output.annot} | cut -f18 | sort -u > {output.promoTargets}
        awk '$8 == "Promoter"' {output.annot} | cut -f1,2,3,6,7 | tail -n +2 > {output.promoBed} # the tail -n +2 is to remove the first line (header) from the output bed
        awk '$8 != "Promoter"' {output.annot} | cut -f1,2,3,6,7 | tail -n +2 > {output.distalBed}
        """

#######################################################################################################################
### BAM TO BIGWIG WITH DEEPTOOLS, GC BIAS
#######################################################################################################################
if NOSPIKE_SAMPLES:
    rule bam2bigwig:
        input: case = lambda wildcards: NOSPIKE_BAM[wildcards.sample], reference = "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/input_{control}.bam"
        output:  "05bigwig/{sample}_{control}-input.bw"
        log: "00log/{sample}_{control}-input_bigwig.bam2bw"
        threads: CLUSTER["bam2bigwig"]["cpu"]
        message: "making input subtracted bigwig for sample {wildcards.sample} with input {input.reference}"
        shell:
            """
            set +u; source activate DPbase; set -u

            inputNorm=$(samtools view -c {input.reference} | awk '{{printf "%f", 1000000*(1/$1)}}') # Scale factor for the input
            sampleNorm=$(samtools view -c {input.case} | awk '{{printf "%f", 1000000*(1/$1)}}') # Scale factor used to scale the sample to substract the input. This step is done to scale with the same method both input and the sample with spike-in to perform the susbtraction.
            printf "%s %f\\n" "scaleFactor for {input.case} is: " "$sampleNorm"

            # Substract the input to the treatment sample. Extend reads, smoothLength and select bin size.
            bamCompare -b1 {input.case} -b2 {input.reference} -o {output} -bs 50 --ratio subtract --scaleFactors $sampleNorm:$inputNorm --extendReads 200
            """

if SPIKE_SAMPLES:
    rule bam2bigwig_spike:
        input: case = lambda wildcards: SPIKE_BAM[wildcards.sample],
                dm = "02aln_dm6/{sample}_dm6.bam",
                reference = "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/input_{control}.bam"
        output: bw = "05bigwig/{sample}_{control}-input.bw", bdg = temp("05bigwig/{sample}_{control}-input.bdg")
        params:
            chr_sizes = "/hpcnfs/scratch/DP/sjammula/mm9/mm9.chrom.sizes",
            bdg2bw = "/hpcnfs/scratch/DP/sjammula/scripts/Tools/bedGraphToBigWig"
        log: "00log/{sample}_{control}-input_bigwig.bam2bw"
        threads: CLUSTER["bam2bigwig"]["cpu"]
        message: "making spike-normalized input subtracted bigwig for sample {wildcards.sample} with input {input.reference}"
        shell:
            """
            set +u; source activate DPbase; set -u

            inputNorm=$(samtools view -c {input.reference} | awk '{{printf "%f", 1000000*(1/$1)}}') # Scale factor for the input
            spikeNorm=$(samtools view -c {input.dm} | awk '{{printf "%f", 1000000*(1/$1)}}') # Scale factor for the spikein sample, used to normalized the final intensity
            sampleNorm=$(samtools view -c {input.case} | awk '{{printf "%f", 1000000*(1/$1)}}') # Scale factor used to scale the sample to substract the input. This step is done to scale with the same method both input and the sample with spike-in to perform the susbtraction.

            # Substract the input to the treatment sample. Extend reads, smoothLength and select bin size.
            bamCompare -b1 {input.case} -b2 {input.reference} -o {output.bw} -bs 50 --ratio subtract --scaleFactors $sampleNorm:$inputNorm --extendReads 200

            # Now we have to revert the previous normalization and apply the spikein factor. So, we multiply the inverse of the sampleNorm by the spikeFactor. This new factor is applied to the bigwig file created.
            sampleNorm2spikeNorm=$(echo $spikeNorm | awk -v sN=$sampleNorm '{{printf "%.2f", $1*(1/sN)}}')
            printf "%s %.2f\\n" "scaleFactor for sampleNorm2spikeNorm is: " "$sampleNorm2spikeNorm"

            # We apply the new normFactor. The output of wiggletools is wig or bedgraph. Then we have to convert again to bigwig.
            wiggletools write_bg {output.bdg} scale $sampleNorm2spikeNorm {output.bw}

            # Transform bedgraph to biwig and finish
            {params.bdg2bw} {output.bdg} {params.chr_sizes} {output.bw}
            """

rule GC_bias:
    input: bam="02aln/{sample}.sorted.bam", bed="03peak_macs2/{sample}_{control}-input/{sample}_peaks_p10.bed"
    output: pdf="06gcBias/{sample}_{control}-input_GCbias.pdf", tmp_txt=temp("06gcBias/{sample}_{control}-input_GCbias.txt")
    log: "00log/{sample}_{control}-input_GCbias.log"
    params: repeatMasker = "/hpcnfs/data/DP/Databases/RepeatMasker_noRandom.bed",
            sumTotalBases = "awk -F\\t 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}'",
            tempBed = temp("06gcBias/{sample}_Repeatmasker.bed")
    message: "Computing GC bias for sample {wildcards.sample}"
    shell:
        """
        set +u; source activate DPbase; set -u
        bedops -u {input.bed} {params.repeatMasker} > {params.tempBed}
        covered=$(bedops --merge {params.tempBed} | {params.sumTotalBases})
        eGs=$((2150570000-$covered))
        computeGCBias -b {input.bam} --effectiveGenomeSize $eGs -g /hpcnfs/data/DP/Databases/mm9.2bit \
        -l 200 -bl {params.tempBed} --biasPlot {output.pdf} --GCbiasFrequenciesFile {output.tmp_txt} 2> {log}
        """

#######################################################################################################################
### FINAL HTML REPORT
#######################################################################################################################
rule multiQC:
    input :
        expand("00log/{sample}.align", sample = ALL_SAMPLES),
        expand("01fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES),
        expand("02aln/{sample}.sorted.bam.flagstat", sample = ALL_SAMPLES)
    output: "10multiQC/multiQC_log.html"
    params: log_name = "multiQC_log"
    log: "00log/multiqc.log"
    message: "multiqc for all logs"
    shell:
        """
        set +u; source activate DPbase; set -u
        multiqc {input} -o 10multiQC -f -v -n {params.log_name} 2> {log}
        """
