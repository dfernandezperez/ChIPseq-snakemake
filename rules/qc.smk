# ------- FASTQC ------- #
rule fastqc:
    input:  
        get_trimmed_forward
    output: 
        "01qc/fqc/{sample}_fastqc.zip"
    log:    
        "00log/fqc/{sample}.log"
    params:
        folder_name = "01qc/fqc/",
        tmp = "01qc/fqc/{sample}.fastq"
    threads: 
        CLUSTER["fastqc"]["cpu"]
    message: 
        "Running fastqc for {input}"
    benchmark:
        ".benchmarks/{sample}.fastqc.benchmark.txt"
    shadow: 
        "minimal"
    shell:
        """
        ln -s "$(pwd)/{input}" {params.tmp}
        fastqc -o {params.folder_name} -f fastq -t {threads} --noextract {params.tmp} 2> {log}
        """


# ------- PhantomPeakQual ------- #
rule phantom_peak_qual:
    input: 
        "02aln/{sample}.bam"
    output:
        "01qc/phantompeakqual/{sample}.spp.out"
    log:
        "00log/phantompeakqual/{sample}_phantompeakqual.log"
    threads:
        CLUSTER["phantom_peak_qual"]["cpu"]
    params:
        out_dir = "01qc/phantompeakqual"
    message:
        "Running phantompeakqual for {wildcards.sample}"
    benchmark:
        ".benchmarks/{sample}.phantompeakqual.benchmark.txt"
    shell:
        """
        /opt/miniconda2/bin/Rscript --vanilla scripts/run_spp_nodups.R \
        -c={input[0]} -savp -rf -p={threads} -odir={params.out_dir}  -out={output} -tmpdir={params.out_dir}  2> {log}
        """


# ------- InsertSize calculation ------- #
rule insert_size:
    input:
        "02aln/{sample}.bam"
    output:
        txt="01qc/insert_size/{sample}.isize.txt",
        pdf="01qc/insert_size/{sample}.isize.pdf"
    log:
        "00log/picard/insert_size/{sample}.log"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE"
    shell:
        """
        # Create the outfiles to handle
        touch {output}
        picard CollectInsertSizeMetrics {params} \
        INPUT={input} OUTPUT={output.txt} \
        HISTOGRAM_FILE={output.pdf} > {log}
        """


# ------- Deeptools quality control ------- #
rule plotFingerprint:
    input: 
        case      = "02aln/{sample}.bam",
        reference = lambda wildcards: "/hpcnfs/data/DP/ChIPseq/INPUT_BAM_FILES/{genome}/input_{control}.bam".format(
            genome = SAMPLES.GENOME[wildcards.sample], 
            control = wildcards.control) 
    output: 
        qualMetrics = "01qc/fingerPrint/{sample}_{control}.qualityMetrics.tsv",
        raw_counts  = "01qc/fingerPrint/{sample}_{control}.rawcounts.tsv",
        plot        = "01qc/fingerPrint/{sample}_{control}.plot.pdf",
    log:
        "00log/plotFingerprint/{sample}_{control}.log"
    params:
        read_exten = set_read_extension,
    threads:
        CLUSTER["plotFingerprint"]["cpu"]
    shell:
        """
        plotFingerprint -b {input} \
        -p {threads} \
        --outQualityMetrics {output.qualMetrics} \
        --outRawCounts {output.raw_counts} \
        --plotFile {output.plot}
        """

rule GC_bias:
    input: 
        bam = "02aln/{sample}.bam",
        bed = rules.filter_peaks.output.bed_filt
    output: 
        pdf      = "01qc/GCbias/{sample}_{control}-input_GCbias.pdf",
        freq_txt = "01qc/GCbias/{sample}_{control}-input_GCbias.txt"
    log:
        "00log/GCbias/{sample}_{control}-input_GCbias.log"
    params:
        repeatMasker = config["ref"]['rep_masker'],
        tempBed      = "01qc/GCbias/{sample}_Repeatmasker.bed.tmp",
        bit_file     = config["ref"]["2bit"],
        egenome_size = config["ref"]["egenome_size"]
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

# ---------------- MultiQC report ----------------- #
rule multiQC_inputs:
    input:
        expand("00log/alignments/{sample}.log", sample = ALL_SAMPLES),
        expand("01qc/fqc/{sample}_fastqc.zip", sample = ALL_SAMPLES),
        expand("01qc/insert_size/{sample}.isize.txt", sample = ALL_SAMPLES),
        expand("01qc/phantompeakqual/{sample}.spp.out", sample = ALL_SAMPLES),
        expand("00log/alignments/rm_dup/{sample}.log", sample = ALL_SAMPLES),
        expand("01qc/fingerPrint/{sample}_{control}.qualityMetrics.tsv", sample = ALL_SAMPLES, control = ALL_CONTROLS),
        expand("01qc/fingerPrint/{sample}_{control}.rawcounts.tsv", sample = ALL_SAMPLES, control = ALL_CONTROLS),
        expand("03peak_macs2/{sample}_{control}-input/{sample}_peaks.xls", sample = ALL_SAMPLES, control = ALL_CONTROLS)
    output: 
        file = "01qc/multiqc_inputs.txt"
    message:
        "create file containing all multiqc input files"
    run:
        with open(output.file, 'w') as outfile:
            for fname in input:
                    outfile.write(fname + "\n")

rule multiQC:
    input:
        "01qc/multiqc_inputs.txt"
    output: 
        "01qc/multiqc_report.html"
    params:
        log_name = "multiqc_report",
        folder = "01qc"
    log:
        "00log/multiqc.log"
    message:
        "multiqc for all logs"
    shell:
        """
        multiqc {input} -o {params.folder} -l {input} -f -v -n {params.log_name} 2> {log}
        """
