# ------- FASTQC ------- #
rule fastqc:
    """This rule is a little bit tricky to fit PE fastq in multiqc.
    Basically the probelm is that fastqc needs to run R1 and R2 separatelly,
    which means 2 fastqc_zip files with different names. This will be recognized
    by multiqc as different samples, so the report will be a mess.
    My workaround has been to use just the forward fastq to create the report.
    For this I need to change the fastq file name (because it has .1.) to fit
    what multiqc is expecting as name. If multiqc reads A.1.fastq it won't know
    that that file must match A.bam in the report, and they will be in different
    rows. I know it's a pain of workaround but it's the only solution I found.
    For this I need to create a symlink to the fastq to be able to change it's name
    (without duplicating the file which would be less efficient). To make sure that 
    the symlink will work it needs to be created from the foldere where it's going to be,
    that's why the cd command of the rule it's imporant. Since the fastq folder can change
    this step needs to work always, it's the only solution I cam up with.
    """
    input:  
        get_trimmed_forward
    output: 
        "results/01qc/fqc/{sample}_fastqc.zip"
    log:    
        "results/00log/fqc/{sample}.log"
    params:
        folder_name = "results/01qc/fqc/",
        tmp = "results/01qc/fqc/{sample}.fastq"
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
        cd {params.folder_name} # Move to folder where symlink is going to be created
        ln -s {input} {params.tmp} # Create symlink to fastq file. Imporant to set the desired file name.
        cd - # Go back to workdir
        fastqc -o {params.folder_name} -f fastq -t {threads} --noextract {params.tmp} 2> {log}
        """


# ------- PhantomPeakQual ------- #
rule phantom_peak_qual:
    input: 
        "results/02aln/{sample}.bam"
    output:
        "results/01qc/phantompeakqual/{sample}.spp.out"
    log:
        "results/00log/phantompeakqual/{sample}_phantompeakqual.log"
    threads:
        CLUSTER["phantom_peak_qual"]["cpu"]
    params:
        out_dir = "results/01qc/phantompeakqual"
    message:
        "Running phantompeakqual for {wildcards.sample}"
    benchmark:
        ".benchmarks/{sample}.phantompeakqual.benchmark.txt"
    shell:
        """
        /opt/miniconda2/bin/Rscript --vanilla workflow/scripts/run_spp_nodups.R \
        -c={input[0]} -savp -rf -p={threads} -odir={params.out_dir}  -out={output} -tmpdir={params.out_dir}  2> {log}
        """


# ------- InsertSize calculation ------- #
rule insert_size:
    input:
        "results/02aln/{sample}.bam"
    output:
        txt="results/01qc/insert_size/{sample}.isize.txt",
        pdf="results/01qc/insert_size/{sample}.isize.pdf"
    log:
        "results/00log/picard/insert_size/{sample}.log"
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
        case      = "results/02aln/{sample}.bam",
        reference = "results/02aln/{control}.bam", 
    output: 
        qualMetrics = "results/01qc/fingerPrint/{sample}_{control}.qualityMetrics.tsv",
        raw_counts  = "results/01qc/fingerPrint/{sample}_{control}.rawcounts.tsv",
        plot        = "results/01qc/fingerPrint/{sample}_{control}.plot.pdf",
    log:
        "results/00log/plotFingerprint/{sample}_{control}.log"
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
        bam = "results/02aln/{sample}.bam",
        bed = rules.filter_peaks.output.bed_filt
    output: 
        pdf      = "results/01qc/GCbias/{sample}_{control}_GCbias.pdf",
        freq_txt = "results/01qc/GCbias/{sample}_{control}_GCbias.txt"
    log:
        "results/00log/GCbias/{sample}_{control}_GCbias.log"
    params:
        repeatMasker = config["ref"]['rep_masker'],
        tempBed      = "results/01qc/GCbias/{sample}_Repeatmasker.bed.tmp",
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
        expand("results/00log/alignments/{sample}.log", sample = SAMPLES.NAME),
        expand("results/01qc/fqc/{sample}_fastqc.zip", sample = SAMPLES.NAME),
        expand("results/01qc/insert_size/{sample}.isize.txt", sample = SAMPLES.NAME),
        expand("results/01qc/phantompeakqual/{sample}.spp.out", sample = ALL_IP),
        expand("results/00log/alignments/rm_dup/{sample}.log", sample = SAMPLES.NAME),
        expand("results/01qc/fingerPrint/{sample}_{control}.qualityMetrics.tsv", zip, sample = ALL_IP, control = ALL_CONTROLS),
        expand("results/01qc/fingerPrint/{sample}_{control}.rawcounts.tsv", zip, sample = ALL_IP, control = ALL_CONTROLS),
        expand("results/03peak_macs2/{sample}_{control}/{sample}_peaks.xls", zip, sample = ALL_IP, control = ALL_CONTROLS)
    output: 
        file = "results/01qc/multiqc/multiqc_inputs.txt"
    message:
        "create file containing all multiqc input files"
    run:
        with open(output.file, 'w') as outfile:
            for fname in input:
                    outfile.write(fname + "\n")

rule multiQC:
    input:
        "results/01qc/multiqc/multiqc_inputs.txt"
    output: 
        "results/01qc/multiqc/multiqc_report.html"
    params:
        log_name = "multiqc_report",
        folder   = "results/01qc/multiqc"
    log:
        "results/00log/multiqc/multiqc.log"
    message:
        "multiqc for all logs"
    shell:
        """
        multiqc -o {params.folder} -l {input} -f -v -n {params.log_name} 2> {log}
        """