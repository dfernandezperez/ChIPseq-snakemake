import pandas as pd
singularity: "/hpcnfs/data/DP/Singularity/dfernandezperez-bioinformatics-singularity-master-chipseq.simg"

#######################################################################################################################
### Load sample sheet and cluster configuration, config file
#######################################################################################################################
configfile: "config.yaml"
CLUSTER = json.load(open(config['CLUSTER_JSON']))
SAMPLES = pd.read_table(config['samples']).set_index("NAME", drop=False)

#######################################################################################################################
### DETERMINE ALL THE OUTPUT FILES TO RUN SNAKEMAKE
#######################################################################################################################
# Get the names of the input files and genome version for calling peak and creating bigwig files
ALL_SAMPLES = SAMPLES.NAME
CONTROLS    = SAMPLES.INPUT
CONTROLS_G  = SAMPLES.GENOME

ALL_CONTROLS  = expand("/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/{genome}/input_{control}.bam", genome = CONTROLS_G, control = CONTROLS)
ALL_FASTQC    = expand("01fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES)
ALL_BAM       = expand("02aln/{sample}.bam", sample = ALL_SAMPLES)
ALL_FLAGSTAT  = expand("02aln/{sample}.bam.flagstat", sample = ALL_SAMPLES)
ALL_PEAKS     = expand("03peak_macs2/{sample}_{control}-input/{sample}_peaks.narrowPeak", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_PHANTOM   = expand("04phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
ALL_PEAKANNOT = expand("05peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + ".annot", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_BIGWIG    = expand( "06bigwig/{sample}_{control}-input.bw", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_GCBIAS    = expand("07gcBias/{sample}_{control}-input_GCbias.pdf", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_QC        = ["10multiQC/multiQC_log.html"]
ALL_BW2SERVER = expand("temp_file_{sample}_{control}.txt",  zip, sample=ALL_SAMPLES, control=CONTROLS)

#######################################################################################################################
### DEFINE LOCAL RULES TO RUN THE WHOLE PIPELINE OR JUST A SUBSET OF IT
#######################################################################################################################
localrules: all, all_server, all_pannot, all_peak_calling, all_fastqc, all_bigwig, all_noGC

rule all:
    input: ALL_FASTQC + ALL_BAM + ALL_FLAGSTAT + ALL_PEAKS + ALL_PHANTOM + ALL_BIGWIG + ALL_QC + ALL_GCBIAS + ALL_PEAKANNOT

rule all_noGC:
    input: ALL_FASTQC + ALL_BAM + ALL_FLAGSTAT + ALL_PEAKS + ALL_PHANTOM + ALL_BIGWIG + ALL_QC + ALL_PEAKANNOT

rule all_server:
    input: ALL_FASTQC + ALL_BAM + ALL_FLAGSTAT + ALL_PEAKS + ALL_PHANTOM + ALL_BIGWIG + ALL_QC + ALL_GCBIAS + ALL_PEAKANNOT + ALL_BW2SERVER

rule all_pannot:
    input: ALL_PEAKANNOT + ALL_PHANTOM

rule all_peak_calling:
    input: ALL_PEAKS

rule all_fastqc:
    input: ALL_FASTQC

rule all_bigwig:
    input: ALL_BIGWIG

#######################################################################################################################
### FASTQC, ALIGNMENT + DEDUP + SABM2BAM SORTED, INDEX BAM
#######################################################################################################################
rule merge_fastqs:
    input: 
        lambda wildcards: SAMPLES.FASTQ[wildcards.sample]
    output: 
        temp("fastq/{sample}.fastq")
    log: 
        "00log/{sample}_unzip"
    threads: 
        CLUSTER["merge_fastqs"]["cpu"]
    params: 
        fastp        = "/hpcnfs/data/DP/software/fastp",
        fastp_params = "-j 00log/{sample}_fastp.json -h 00log/{sample}_fastp.html --stdin -t 1 -A -Q -L "
    message: 
        "merging fastq files {input} into {output}"
    benchmark:
        ".benchmarks/{sample}.merge_fastqs.benchmark.txt"
    shell:
        """
        file={input}
        if [ ${{file: -8}} == "fastq.gz" ]
        then
            fastq={input}
        else
            fastq={input}/*fastq.gz
        fi
        zcat $fastq | {params.fastp} -o {output} -w {threads} {params.fastp_params} 2> {log}
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
    benchmark:
        ".benchmarks/{sample}.fastqc.benchmark.txt"
    shell:
        """
        fastqc -o 01fqc -f fastq -t {threads} --noextract {input} 2> {log}
        """

# Check if there's spike in, get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and duplicated reads by samblaster -r
# samblaster should run before samtools sort
rule align:
    input:
        lambda wildcards: 
            str("fastq/" + wildcards.sample + ".fastq") if SAMPLES.SPIKE[wildcards.sample] == False else str()
    output:
        "02aln/{sample}.bam"
    threads:
        CLUSTER["align"]["cpu"]
    params:
        bowtie = "--chunkmbs 1024 -m 1 --best -S " + config["idx_bt1_mm"]
    message:
        "Aligning {input} with parameters {params.bowtie}"
    log:
        bowtie     = "00log/{sample}.align",
        markdup    = "00log/{sample}.markdup",
        sort_index = "00log/{sample}.sortIndex_bam"
    benchmark:
        ".benchmarks/{sample}.align.benchmark.txt"
    shell:
        """
        bowtie -p {threads} {params.bowtie} -q {input} 2> {log.bowtie} \
        | samblaster --removeDups 2> {log.markdup} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m 2G -@ {threads} -T {output}.tmp -o {output} - 2> {log.sort_index}
        samtools index {output} 2>> {log.sort_index}
        """

rule align_spike:
    input:
        lambda wildcards: 
            str("fastq/" + wildcards.sample + ".fastq") if SAMPLES.SPIKE[wildcards.sample] == True else str()
    output:
        mm = "02aln/{sample}.bam",
        dm = "02aln_dm/{sample}_dm.bam"
    threads:
        CLUSTER["align"]["cpu"]
    params:
        bowtie_mm = "--chunkmbs 1024 -m 1 --best -S " + config["idx_bt1_mm"],
        bowtie_dm = "--chunkmbs 1024 -m 1 --best -S " + config["idx_bt1_dm"]
    message:
        "Aligning {input} with parameters {params.bowtie_mm}"
    log:
        bowtie     = ["00log/{sample}.align", "00log/{sample}.dm_align"],
        markdup    = "00log/{sample}.markdup",
        sort_index = "00log/{sample}.sortIndex_bam"
    benchmark:
        ".benchmarks/{sample}.align-spike.benchmark.txt"
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
        | samtools sort -m 2G -@ {threads} -T {output.dm}.tmp -o {output.dm} - 2>> {log.sort_index}
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
    benchmark:
        ".benchmarks/{sample}.flagstat.benchmark.txt"
    shell:
        """
        samtools flagstat {input} > {output} 2> {log}
        """

#######################################################################################################################
### PEAK CALLING WITH MACS2 AND PEAK QUALITY STATS WITH SPP
#######################################################################################################################
rule call_peaks:
    input: 
        case      = "02aln/{sample}.bam",
        reference = lambda wildcards: "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/" + SAMPLES.GENOME[wildcards.sample] + "/input_" + wildcards.control + ".bam",
    output: 
        narrowPeak = "03peak_macs2/{sample}_{control}-input/{sample}_peaks.narrowPeak"
    log:
        "00log/{sample}_{control}-input_macs2.log"
    params:
        out_dir      = "03peak_macs2/{sample}_{control}-input",
        macs2_params = "--keep-dup all --format BAM --nomodel -m 10 30",
        pvalue       = config["macs2"]["pvalue"],
        gsize        = config["macs2"]["gsize"]
    message: 
        "call_peaks macs2 with input {input.reference} for sample {input.case}"
    benchmark:
        ".benchmarks/{sample}_{control}.callpeaks.benchmark.txt"
    shell:
        """
        macs2 callpeak --treatment {input.case} \
            --control {input.reference} \
            --gsize {params.gsize} \
            --outdir {params.out_dir} \
            --name {wildcards.sample} \
            --pvalue {params.pvalue} \
            {params.macs2_params} 2> {log}            
        """

rule filter_peaks:
    input:
        rules.call_peaks.output.narrowPeak
    output:
        bed_filt = "03peak_macs2/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + ".bed"
    params:
        pval_filt = config["macs2"]["filt_peaks_pval"]
    shell:
        """
        awk "\$8 >= {params.pval_filt}" {input} | cut -f1-4,8 > {output.bed_filt}
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
    benchmark:
        ".benchmarks/{sample}.phantompeakqual.benchmark.txt"
    shell:
        """
        /opt/miniconda2/bin/Rscript --vanilla scripts/run_spp_nodups.R \
        -c={input[0]} -savp -rf -p={threads} -odir={params.out_dir}  -out={output} -tmpdir={params.out_dir}  2> {log}
        """

rule peakAnnot:
    input :
        rules.filter_peaks.output.bed_filt,
    output:
        annot             = "05peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + ".annot",
        promo_bed_targets = "05peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + "_promoTargets.bed",
        promoTargets      = "05peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + "_promoTargets.txt",
        promoBed          = "05peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + "_promoPeaks.bed",
        distalBed         = "05peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + "_distalPeaks.bed"
    params:
        before = config["promoter"]["bTSS"],
        after  = config["promoter"]["aTSS"],
        genome = lambda wildcards: SAMPLES.GENOME[wildcards.sample]
    log: 
        "00log/{sample}_{control}-input_peakanot"
    message:
        "Annotating peaks for {wildcards.sample}"
    shell:
        """
        Rscript --vanilla scripts/PeakAnnot.R {input} {params.before} {params.after}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} {output.distalBed} {params.genome}
        """

#######################################################################################################################
### BAM TO BIGWIG WITH DEEPTOOLS, BDG TO BIGWIG TO SERVER, GC BIAS
#######################################################################################################################
rule bam2bigwig:
    input: 
        case      = lambda wildcards: str("02aln/" + wildcards.sample + ".bam") if SAMPLES.SPIKE[wildcards.sample] == False else str(),
        reference = lambda wildcards: "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/" + SAMPLES.GENOME[wildcards.sample] + "/input_" + wildcards.control + ".bam",
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
        case      = lambda wildcards: str("02aln/" + wildcards.sample + ".bam") if SAMPLES.SPIKE[wildcards.sample] == True else str(),
        reference = lambda wildcards: "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/" + SAMPLES.GENOME[wildcards.sample] + "/input_" + wildcards.control + ".bam",
        dm        = "02aln_dm/{sample}_dm.bam"
    output:
        "06bigwig/{sample}_{control}-input.bw"
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

rule bigwig2server:
    input: 
        bw       = "06bigwig/{sample}_{control}-input.bw",
        flagstat = "02aln/{sample}.bam.flagstat"
    output:
        temp("temp_file_{sample}_{control}.txt")
    params:
        user     = lambda wildcards : SAMPLES.USER[wildcards.sample],
        antibody = lambda wildcards : SAMPLES.AB[wildcards.sample],
        genome   = lambda wildcards : SAMPLES.GENOME[wildcards.sample],
        run      = lambda wildcards : SAMPLES.RUN[wildcards.sample],
        chip     = lambda wildcards : str("ChIPseq") if SAMPLES.SPIKE[wildcards.sample] == False else str("ChIPseqSpike")
    run:
        with open (input.flagstat, "r") as f:
            line         = f.readlines()[4] # fifth line contains the number of mapped reads
            match_number = re.match(r'(\d.+) \+.+', line)
            total_reads  = int(match_number.group(1))
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

rule GC_bias:
    input: 
        bam = "02aln/{sample}.bam",
        bed = rules.filter_peaks.output.bed_filt
    output: 
        pdf      = "07gcBias/{sample}_{control}-input_GCbias.pdf",
        freq_txt = "07gcBias/{sample}_{control}-input_GCbias.txt"
    log:
        "00log/{sample}_{control}-input_GCbias.log"
    params:
        repeatMasker = config["genome"]['rep_masker'],
        tempBed      = "07gcBias/{sample}_Repeatmasker.bed.tmp",
        bit_file     = config["genome"]["2bit"],
        egenome_size = config["genome"]["egenome_size"]
    threads:
        CLUSTER["GC_bias"]["cpu"]
    message:
        "Computing GC bias for sample {wildcards.sample}"
    benchmark:
        ".benchmarks/{sample}_{control}.GCbias.benchmark.txt"
    shell:
        """
        bedops -u {input.bed} {params.repeatMasker} > {params.tempBed}
        bp_peaks=$(bedops --merge {input.bed} | bedmap --bases - | awk "{{sum+=\$1}}END{{print sum}}")
        total_eGsize=$(({params.egenome_size}-$bp_peaks))

        computeGCBias -b {input.bam} \
            -p {threads} \
            --effectiveGenomeSize $total_eGsize \
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