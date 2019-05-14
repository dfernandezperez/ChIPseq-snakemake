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
ALL_PHANTOM   = expand("04phantompeakqual/{sample}_phantom.txt", sample = ALL_SAMPLES)
ALL_PEAKANNOT = expand("05peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["macs2"]["filt_peaks_pval"] + ".annot", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_BIGWIG    = expand( "06bigwig/{sample}_{control}-input.bw", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_GCBIAS    = expand("07gcBias/{sample}_{control}-input_GCbias.pdf", zip, sample=ALL_SAMPLES, control=CONTROLS)
ALL_QC        = ["10multiQC/multiQC_log.html"]
# ALL_BW2SERVER = expand("temp_file_{sample}_{control}.txt",  zip, sample=ALL_SAMPLES, control=CONTROLS)

#######################################################################################################################
### DEFINE LOCAL RULES TO RUN THE WHOLE PIPELINE OR JUST A SUBSET OF IT
#######################################################################################################################
localrules: all, all_server, all_pannot, all_peak_calling, all_fastqc, all_bigwig, all_noGC

rule all:
    input: ALL_FASTQC + ALL_BAM + ALL_FLAGSTAT + ALL_PEAKS + ALL_PHANTOM + ALL_BIGWIG + ALL_QC + ALL_GCBIAS + ALL_PEAKANNOT

rule all_noGC:
    input: ALL_FASTQC + ALL_BAM + ALL_FLAGSTAT + ALL_PEAKS + ALL_PHANTOM + ALL_BIGWIG + ALL_QC + ALL_PEAKANNOT

rule all_server:
    input: ALL_FASTQC + ALL_BAM + ALL_FLAGSTAT + ALL_PEAKS + ALL_PHANTOM + ALL_BIGWIG + ALL_QC + ALL_PEAKANNOT + ALL_BW2SERVER

rule all_pannot:
    input: ALL_PEAKANNOT + ALL_PHANTOM

rule all_peak_calling:
    input: ALL_PEAKS

rule all_fastqc:
    input: ALL_FASTQC

rule all_bigwig:
    input: ALL_BIGWIG



##### load rules #####

include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/peaks.smk"
include: "rules/qc.smk"


##### handle possible errors, clean temp folders #####
onsuccess:
    shell("rm -r fastq/")

onerror:
    print("An error ocurred. Workflow aborted")
    shell("""
        mail -s "An error occurred. ChIP-seq snakemake workflow aborted" `whoami`@ieo.it < {log}
        """)
