rule phantom_peak_qual:
    input: 
        get_bam
    output:
        "04phantompeakqual/{sample}_phantom.txt"
    log:
        "00log/qc/{sample}_phantompeakqual.log"
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


rule GC_bias:
    input: 
        bam = get_bam,
        bed = rules.filter_peaks.output.bed_filt
    output: 
        pdf      = "07gcBias/{sample}_{control}-input_GCbias.pdf",
        freq_txt = "07gcBias/{sample}_{control}-input_GCbias.txt"
    log:
        "00log/qc/{sample}_{control}-input_GCbias.log"
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

# ---------------- MultiQC report ----------------- #
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