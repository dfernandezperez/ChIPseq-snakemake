rule call_peaks:
    input: 
        case      = "results/02aln/{sample}.bam",
        reference = "results/02aln/{control}.bam"
    output: 
        narrowPeak = "results/03peak_macs2/{sample}_{control}/{sample}_peaks.narrowPeak",
	summit     = "results/03peak_macs2/{sample}_{control}/{sample}_summits.bed",
        xls        = "results/03peak_macs2/{sample}_{control}/{sample}_peaks.xls"
    log:
        "results/00log/macs2/{sample}_{control}_macs2.log"
    params:
        out_dir      = "results/03peak_macs2/{sample}_{control}",
        macs2_params = config["params"]["macs2"]["pk_calling"],
        pvalue       = config["params"]["macs2"]["pvalue"],
        gsize        = config["params"]["macs2"]["gsize"],
        paired_end = lambda w: "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
    message: 
        "call_peaks macs2 with input {input.reference} for sample {input.case}"
    benchmark:
        "results/.benchmarks/{sample}_{control}.callpeaks.benchmark.txt"
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


rule call_peaks_broad:
    input: 
        case      = "results/02aln/{sample}.bam",
        reference = "results/02aln/{control}.bam"
    output: 
        broadPeak  = "results/03peak_macs2/{sample}_{control}/broad/{sample}_peaks.broadPeak",
        gappedPeak = "results/03peak_macs2/{sample}_{control}/broad/{sample}_peaks.gappedPeak",
        xls        = "results/03peak_macs2/{sample}_{control}/broad/{sample}_peaks.xls",
    log:
        "results/00log/macs2/{sample}_{control}_macs2.log"
    params:
        out_dir       = "results/03peak_macs2/{sample}_{control}/broad",
        macs2_params  = config["params"]["macs2"]["pk_calling"],
        pvalue        = config["params"]["macs2"]["pvalue_broad"],
        pvalue_narrow = config["params"]["macs2"]["pvalue"],
        gsize         = config["params"]["macs2"]["gsize"],
        paired_end    = lambda w: "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
    message: 
        "call_peaks macs2 with input {input.reference} for sample {input.case}"
    benchmark:
        "results/.benchmarks/{sample}_{control}.callpeaks.benchmark.txt"
    shell:
        """
        macs2 callpeak {params.paired_end} \
            --broad --broad-cutoff {params.pvalue} \
            --treatment {input.case} \
            --control {input.reference} \
            --gsize {params.gsize} \
            --outdir {params.out_dir} \
            --name {wildcards.sample} \
            --pvalue {params.pvalue_narrow} \
            {params.macs2_params} 2> {log}            
        """


rule filter_peaks:
    input:
        narrowpeak = rules.call_peaks.output.narrowPeak
	summit     = rules.call_peaks.output.summit
    output:
        bed_filt    = "results/03peak_macs2/{{sample}}_{{control}}/{{sample}}_peaks_p{pvalue}.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
	summit_filt = "results/03peak_macs2/{{sample}}_{{control}}/{{sample}}_summits_p{pvalue}.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
    params:
        pval_filt = config["params"]["macs2"]["filt_peaks_pval"]
    shell:
        """
        awk "\$8 >= {params.pval_filt}" {input.narrowpeak} | cut -f1-4,8 > {output.bed_filt}
        awk "\$5 >= {params.pval_filt}" {input.summit} > {output.summit_filt}
	"""


rule peakAnnot:
    input :
        rules.filter_peaks.output.bed_filt,
    output:
        annot             = "results/04peak_annot/{{sample}}_{{control}}/{{sample}}_peaks_p{pvalue}.annot".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promo_bed_targets = "results/04peak_annot/{{sample}}_{{control}}/{{sample}}_peaks_p{pvalue}_promoTargets.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoTargets      = "results/04peak_annot/{{sample}}_{{control}}/{{sample}}_peaks_p{pvalue}_promoTargets.txt".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoBed          = "results/04peak_annot/{{sample}}_{{control}}/{{sample}}_peaks_p{pvalue}_promoPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        distalBed         = "results/04peak_annot/{{sample}}_{{control}}/{{sample}}_peaks_p{pvalue}_distalPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
    params:
        before = config["promoter"]["bTSS"],
        after  = config["promoter"]["aTSS"],
        genome = lambda wildcards: SAMPLES.GENOME[wildcards.sample]
    log: 
        "results/00log/peakAnnot/{sample}_{control}_peakanot"
    message:
        "Annotating peaks for {wildcards.sample}"
    shell:
        """
        Rscript --vanilla workflow/scripts/PeakAnnot.R {input} {params.before} {params.after}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} {output.distalBed} {params.genome}
        """
