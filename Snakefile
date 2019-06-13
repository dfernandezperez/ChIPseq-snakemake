import pandas as pd
import yaml
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.4.3")

# Singularity image path
singularity: "/hpcnfs/data/DP/Singularity/dfernandezperez-bioinformatics-singularity-master-chipseq.simg"

#######################################################################################################################
### Load sample sheet and cluster configuration, config file
#######################################################################################################################
configfile: "config.yaml"

CLUSTER     = yaml.load(open(config['cluster'], 'r'))
SAMPLES     = pd.read_table(config['samples']).set_index("NAME", drop=False)
units       = pd.read_table(config["units"], dtype=str).set_index(["sample", "lane"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index


#######################################################################################################################
### DETERMINE ALL THE OUTPUT FILES TO RUN SNAKEMAKE
#######################################################################################################################
# Get the names of the input files and genome version for calling peak and creating bigwig files
ALL_SAMPLES  = SAMPLES.NAME
ALL_CONTROLS = SAMPLES.INPUT
CONTROLS_G   = SAMPLES.GENOME

ALL_PEAKANNOT = expand("04peak_annot/{sample}_{control}-input/{sample}_peaks_p" + config["params"]["macs2"]["filt_peaks_pval"] + ".annot", zip, sample=ALL_SAMPLES, control=ALL_CONTROLS)
ALL_PEAKS     = expand("03peak_macs2/{sample}_{control}-input/{sample}_peaks.narrowPeak", zip, sample=ALL_SAMPLES, control=ALL_CONTROLS)
ALL_BIGWIG    = expand("06bigwig/{sample}_{control}-input.bw", zip, sample=ALL_SAMPLES, control=ALL_CONTROLS)
ALL_GCBIAS    = expand("01qc/GCbias/{sample}_{control}-input_GCbias.pdf", zip, sample=ALL_SAMPLES, control=ALL_CONTROLS)
ALL_QC        = ["01qc/multiqc/multiqc_report.html"]
# ALL_BW2SERVER = expand("temp_file_{sample}_{control}.txt",  zip, sample=ALL_SAMPLES, control=ALL_CONTROLS)


#######################################################################################################################
### DEFINE LOCAL RULES TO RUN THE WHOLE PIPELINE OR JUST A SUBSET OF IT
#######################################################################################################################
localrules: all, all_noGC, test

rule all:
    input: ALL_PEAKANNOT + ALL_BIGWIG + ALL_QC + ALL_GCBIAS

rule all_noGC:
    input:  ALL_PEAKANNOT + ALL_BIGWIG + ALL_QC

# rule all_server:
#     input: ALL_PHANTOM + ALL_PEAKANNOT + ALL_BIGWIG + ALL_QC + ALL_GCBIAS + ALL_BW2SERVER


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