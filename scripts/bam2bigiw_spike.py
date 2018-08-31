import subprocess
import pysam

####################
## Read bam files ##
####################
reference = pysam.AlignmentFile(snakemake.input["reference"], "rb")
dm = pysam.AlignmentFile(snakemake.input["dm"], "rb")
case = pysam.AlignmentFile(snakemake.input["case"], "rb")

################################################################
## Calculate normalization factors: (1/mapped reads)*1million ##
################################################################
reference_norm = (1.0/float(reference.mapped))*1000000
dm_norm = (1.0/float(dm.mapped))*1000000
case_norm = (1.0/float(case.mapped))*1000000

sampleNorm2spikeNorm = (1.0/case_norm)*dm_norm

#############################
## Bash commands to launch ##
#############################
bamCompare = "bamCompare -b1 " + snakemake.input["case"] + " -b2 " + snakemake.input["reference"] + " -o " + snakemake.output["bw"] + " -bs 50 --ratio subtract --scaleFactors " + case_norm + ":" + dm_norm + " --extendReads 200"
wiggleTools = "wiggletools write_bg " + snakemake.output["bdg"] + " scale " + $sampleNorm2spikeNorm + " " + {output.bw}
bdg2bw = "bedGraphToBigWig " + snakemake.output["bdg"] + " " + snakemake.params["chr_sizes"] + " " + snakemake.output["bw"]

subprocess.Popen(bamCompare.split())
subprocess.Popen(wiggleTools.split())
subprocess.Popen(bdg2bw.split())
