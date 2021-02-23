rule set_downsample_proportion:
    input:
        lambda w: expand("results/02aln/{sample}.bam", sample = IPS_downsample[ IPS_downsample["DOWNSAMPLE_GROUP"] == w.group ].NAME)
    output:
        downsample_factor = "results/02aln/downsampled/set_proportion_{group}.txt"
    log:
        "results/00log/set_downsample_proportion/{group}.log"
    shell:
        """
        python workflow/scripts/get_proportion_downsample.py \
        --input {input} \
        --output {output}
        """

rule perform_downsample:
    input:
        bam               = "results/02aln/{sample}.bam",
        downsample_factor = "results/02aln/downsampled/set_proportion_{group}.txt"
    output:
        bam = temp("results/02aln/downsampled/{sample}_{group}.bam")
    shell:
        """
        target=`cat {input.downsample_factor}`
        nreads=`samtools view -c {input.bam}`
        prop=$(echo "print($target/$nreads)" | python3)

        samtools view -bs {input.bam} > {output}
        """

def select_input(wildcards):
    if SAMPLES.loc[ SAMPLES.loc[wildcards.sample].INPUT ].DOWNSAMPLE_GROUP == "FALSE":
        return { "case" : "results/02aln/downsampled/{sample}_{group}.bam",
                "reference" : "results/02aln/{control}.bam" }
    else:
        dwn_group = SAMPLES.loc[ SAMPLES.loc[wildcards.sample].INPUT ].DOWNSAMPLE_GROUP
        return { "case" : "results/02aln/downsampled/{sample}_{group}.bam",
                "reference" : "results/02aln/downsampled/{control}_{group_input}.bam".format(group_input = dwn_group, control = wildcards.control) }


rule call_peaks_downsample:
    input: 
        unpack(select_input)
    output: 
        narrowPeak = "results/03peak_macs2/downsampled/{sample}_{control}_{group}/{sample}_peaks.narrowPeak",
        xls        = "results/03peak_macs2/downsampled/{sample}_{control}_{group}/{sample}_peaks.xls"
    log:
        "results/00log/macs2/downsampled/{sample}_{control}_{group}_macs2.log"
    params:
        out_dir      = "results/03peak_macs2/downsampled/{sample}_{control}_{group}",
        macs2_params = config["params"]["macs2"]["pk_calling"],
        pvalue       = config["params"]["macs2"]["pvalue"],
        gsize        = config["params"]["macs2"]["gsize"],
        paired_end = lambda w: "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
    message: 
        "call_peaks macs2 with input {input.reference} for sample {input.case}"
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


rule call_peaks_broad_downsample:
    input: 
        case      = "results/02aln/downsample/{sample}_{group}.bam",
        reference = "results/02aln/{control}.bam"
    output: 
        broadPeak  = "results/03peak_macs2/downsampled/{sample}_{control}_{group}/broad/{sample}_peaks.broadPeak",
        gappedPeak = "results/03peak_macs2/downsampled/{sample}_{control}_{group}/broad/{sample}_peaks.gappedPeak",
        xls        = "results/03peak_macs2/downsampled/{sample}_{control}_{group}/broad/{sample}_peaks.xls",
    log:
        "results/00log/macs2/downsampled/{sample}_{control}_{group}_macs2_broad.log"
    params:
        out_dir       = "results/03peak_macs2/downsampled/{sample}_{control}_{group}/broad",
        macs2_params  = config["params"]["macs2"]["pk_calling"],
        pvalue        = config["params"]["macs2"]["pvalue_broad"],
        pvalue_narrow = config["params"]["macs2"]["pvalue"],
        gsize         = config["params"]["macs2"]["gsize"],
        paired_end    = lambda w: "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
    message: 
        "call_peaks macs2 with input {input.reference} for sample {input.case}"
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


rule filter_peaks_downsample:
    input:
        rules.call_peaks_downsample.output.narrowPeak
    output:
        bed_filt = "results/03peak_macs2/downsampled/{{sample}}_{{control}}_{{group}}/{{sample}}_peaks_p{pvalue}.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
    params:
        pval_filt = config["params"]["macs2"]["filt_peaks_pval"]
    shell:
        """
        awk "\$8 >= {params.pval_filt}" {input} | cut -f1-4,8 > {output.bed_filt}
        """


rule peakAnnot_downsample:
    input :
        rules.filter_peaks_downsample.output.bed_filt,
    output:
        annot             = "results/04peak_annot/downsampled/{{sample}}_{{control}}_{{group}}/{{sample}}_peaks_p{pvalue}.annot".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promo_bed_targets = "results/04peak_annot/downsampled/{{sample}}_{{control}}_{{group}}/{{sample}}_peaks_p{pvalue}_promoTargets.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoTargets      = "results/04peak_annot/downsampled/{{sample}}_{{control}}_{{group}}/{{sample}}_peaks_p{pvalue}_promoTargets.txt".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        promoBed          = "results/04peak_annot/downsampled/{{sample}}_{{control}}_{{group}}/{{sample}}_peaks_p{pvalue}_promoPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"]),
        distalBed         = "results/04peak_annot/downsampled/{{sample}}_{{control}}_{{group}}/{{sample}}_peaks_p{pvalue}_distalPeaks.bed".format(pvalue = config["params"]["macs2"]["filt_peaks_pval"])
    params:
        before = config["promoter"]["bTSS"],
        after  = config["promoter"]["aTSS"],
        genome = lambda wildcards: SAMPLES.GENOME[wildcards.sample]
    log: 
        "results/00log/peakAnnot/downsampled/{sample}_{control}_{group}_peakanot"
    message:
        "Annotating peaks for {wildcards.sample}"
    shell:
        """
        Rscript --vanilla workflow/scripts/PeakAnnot.R {input} {params.before} {params.after}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} {output.distalBed} {params.genome}
        """