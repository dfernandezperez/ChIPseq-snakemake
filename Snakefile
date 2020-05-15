import pandas as pd
import numpy as np
import yaml
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("5.4.3")

##### Singularity image path
singularity: "/hpcnfs/data/DP/Singularity/dfernandezperez-bioinformatics-singularity-master-chipseq.simg"

##### config file
configfile: "configuration/config.yaml"


#######################################################################################################################
### Load sample sheet and cluster configuration, config file
#######################################################################################################################
SAMPLES     = pd.read_csv(config['samples'], sep ="\t").set_index("NAME", drop=False).sort_index()
units       = pd.read_csv(config["units"], dtype=str, sep ="\t").set_index(["sample", "lane"], drop=False).sort_index()
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index

CLUSTER     = yaml.load(open(config['cluster'], 'r'))


#######################################################################################################################
### DETERMINE ALL THE OUTPUT FILES TO RUN SNAKEMAKE
#######################################################################################################################
# Subset the samples that are IP (no inputs)
# Get the names of the input files and genome version for calling peak and creating bigwig files
IPS          = SAMPLES[ SAMPLES["IS_INPUT"] == False ]
ALL_IP       = IPS.NAME
ALL_CONTROLS = IPS.INPUT
CONTROLS_G   = IPS.GENOME


# Define the output files for the rule ALL
ALL_PEAKANNOT    = expand("04peak_annot/{sample}_{control}/{sample}_peaks_p" + config["params"]["macs2"]["filt_peaks_pval"] + ".annot", zip, sample = ALL_IP, control = ALL_CONTROLS)
ALL_PEAKS        = expand("03peak_macs2/{sample}_{control}/{sample}_peaks.narrowPeak", zip, sample = ALL_IP, control = ALL_CONTROLS)
ALL_BROAD_PEAKS  = expand("03peak_macs2/{sample}_{control}/broad/{sample}_peaks.broadPeak", zip, sample=ALL_IP, control=ALL_CONTROLS)
ALL_BIGWIG       = expand("06bigwig/{sample}_{control}.bw", zip, sample = ALL_IP, control = ALL_CONTROLS)
ALL_GCBIAS       = expand("01qc/GCbias/{sample}_{control}_GCbias.pdf", zip, sample = ALL_IP, control = ALL_CONTROLS)
ALL_QC           = ["01qc/multiqc/multiqc_report.html"]
ALL_BIGWIG_NOSUB = expand("06bigwig/noSubtract/{sample}.bw", sample = SAMPLES.NAME)
ALL_BW2SERVER    = expand("temp_file_{sample}_{control}.txt",  zip, sample = ALL_IP, control = ALL_CONTROLS)


#-------------------- Set variables and target files to prepare data for GEO upload -----------------------#
ALL_IP_SE = set(units[units['fq2'].isnull()]['sample'])
ALL_IP_PE = set(units[units['fq2'].notnull()]['sample'])

ALL_FASTQ_GEO_PE = expand(["GEO/fastq/{sample}.1.fastq.gz", "GEO/fastq/{sample}.2.fastq.gz"], sample = ALL_IP_PE)
ALL_FASTQ_GEO_SE = expand("GEO/fastq/{sample}.se.fastq.gz", sample = ALL_IP_SE)
ALL_PEAKS_GEO    = expand(
    "GEO/peaks/{sample}_{control}_peaks_p" + config["params"]["macs2"]["filt_peaks_pval"] + ".bed", 
    zip, 
    sample = ALL_IP,
    control = ALL_CONTROLS
    )


###########################################################################################################################
### Local rules are rules that won't be submitted to the scheduler but executed in the current session (front-end or node)
###########################################################################################################################
localrules: cp_fastq_pe, cp_fastq_se, all, all_broad, GC, server, geo


###########################################################################################################################
### Define multiple outputs based on the output files desired
###########################################################################################################################
rule all:
    input:  ALL_PEAKS + ALL_PEAKANNOT + ALL_BIGWIG + ALL_QC + ALL_BIGWIG_NOSUB

rule all_broad:
    input:  ALL_PEAKS + ALL_PEAKANNOT + ALL_BIGWIG + ALL_QC + ALL_BIGWIG_NOSUB + ALL_BROAD_PEAKS

rule GC:
    input:  ALL_GCBIAS

rule server:
	input:  ALL_BW2SERVER

rule geo:
    input: "GEO/md5sum/md5sum_peaks.txt", "GEO/md5sum/md5sum_fastqs.txt", ALL_PEAKS_GEO + ALL_FASTQ_GEO_SE + ALL_FASTQ_GEO_PE


##### load rules #####
include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/peaks.smk"
include: "rules/qc.smk"
include: "rules/prepare2GEO.smk"

##### handle possible errors, clean temp folders #####
# Remove the folder used to create the fastq files (snakemake removes the tmp files but not the folder...)
# Since some jobs a lot of times end in E state after finishing (when they're too fast, like creating a soft link),
# remove those "canceled" jobs after the pipeline ends
onsuccess:
    shell(
	"""
    rm -r fastq/
    qselect -u `whoami` -s E | xargs qdel -Wforce
    """)

onerror:
    print("An error ocurred. Workflow aborted")
    shell(
	"""
	qselect -u `whoami` -s E | xargs qdel -Wforce
    mail -s "An error occurred. ChIP-seq snakemake workflow aborted" `whoami`@ieo.it < {log}
    """)
