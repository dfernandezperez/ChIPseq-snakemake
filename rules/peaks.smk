rule call_peaks:
    input: 
        case      = "02aln/{sample}.bam",
        reference = lambda wildcards: "/hpcnfs/scratch/DP/dfernand/Stamburri/PCGFS_PROJECT/ChIPseq/INPUT/{genome}/input_{control}.bam".format(
            genome = SAMPLES.GENOME[wildcards.sample], 
            control = wildcards.control) 
    output: 
        narrowPeak = "03peak_macs2/{sample}_{control}-input/{sample}_peaks.narrowPeak",
        xls        = "03peak_macs2/{sample}_{control}-input/{sample}_peaks.xls"
    log:
        "00log/macs2/{sample}_{control}-input_macs2.log"
    params:
        out_dir      = "03peak_macs2/{sample}_{control}-input",
        macs2_params = config["params"]["macs2"]["pk_calling"],
        pvalue       = config["params"]["macs2"]["pvalue"],
        gsize        = config["params"]["macs2"]["gsize"],
        paired_end = lambda w: "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
    message: 
        "call_peaks macs2 with input {input.reference} for sample {input.case}"
    benchmark:
        ".benchmarks/{sample}_{control}.callpeaks.benchmark.txt"
    shell:
        """
        macs2 callpeak {params.paired_end} \
            --treatment {input.case} \
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
        bed_filt = "03peak_macs2/{{sample}}_{{control}}-input/{{sample}}_peaks_p{pvalue}.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
    params:
        pval_filt = config["params"]["macs2"]["filt_peaks_pval"]
    shell:
        """
        awk "\$8 >= {params.pval_filt}" {input} | cut -f1-4,8 > {output.bed_filt}
        """


rule peakAnnot:
    input :
        rules.filter_peaks.output.bed_filt,
    output:
        annot             = "04peak_annot/{{sample}}_{{control}}-input/{{sample}}_peaks_p{pvalue}.annot".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promo_bed_targets = "04peak_annot/{{sample}}_{{control}}-input/{{sample}}_peaks_p{pvalue}_promoTargets.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoTargets      = "04peak_annot/{{sample}}_{{control}}-input/{{sample}}_peaks_p{pvalue}_promoTargets.txt".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoBed          = "04peak_annot/{{sample}}_{{control}}-input/{{sample}}_peaks_p{pvalue}_promoPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        distalBed         = "04peak_annot/{{sample}}_{{control}}-input/{{sample}}_peaks_p{pvalue}_distalPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
    params:
        before = config["promoter"]["bTSS"],
        after  = config["promoter"]["aTSS"],
        genome = lambda wildcards: SAMPLES.GENOME[wildcards.sample]
    log: 
        "00log/peakAnnot/{sample}_{control}-input_peakanot"
    message:
        "Annotating peaks for {wildcards.sample}"
    shell:
        """
        Rscript --vanilla scripts/PeakAnnot.R {input} {params.before} {params.after}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} {output.distalBed} {params.genome}
        """