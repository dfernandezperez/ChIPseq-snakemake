import subprocess
import pysam

def touch_file(file):
	command = "touch " + file
	subprocess.call(command.split())

# Avoid the warning from pysam that tells that the bam file is more recent than the index file (due to snakemake behaviour)
touch_file(snakemake.input["dm"] + ".bai")
touch_file(snakemake.input["case"] + ".bai")

####################
## Read bam files ##
####################
reference = pysam.AlignmentFile(snakemake.input["reference"], "rb")
dm = pysam.AlignmentFile(snakemake.input["dm"], "rb")
case = pysam.AlignmentFile(snakemake.input["case"], "rb")

################################################################
## Calculate normalization factors: (1/mapped reads)*1million ##
################################################################
reference_norm = str( (1.0/float(reference.mapped))*1000000 )
dm_norm = str( (1.0/float(dm.mapped))*1000000 )
case_norm = str( (1.0/float(case.mapped))*1000000 )

sampleNorm2spikeNorm = str( (1.0/float(case_norm))*float(dm_norm) )

#############################
## Bash commands to launch ##
#############################
bamCompare = "bamCompare -b1 " + snakemake.input["case"] + " -b2 " + snakemake.input["reference"] + " -o " + snakemake.output["bw"] + " -bs 50 --ratio subtract --scaleFactors " + case_norm + ":" + reference_norm + " " + snakemake.threads + " --extendReads 200 &> " + snakemake.log

wiggleTools = "wiggletools write_bg " + snakemake.output["bdg"] + " scale " + sampleNorm2spikeNorm + " " + snakemake.output["bw"]

#bedGraphToBigWig = "/hpcnfs/scratch/DP/sjammula/scripts/Tools/bedGraphToBigWig"
bdg2bw = "bedGraphToBigWig " + snakemake.output["bdg"] + " " + snakemake.params["chr_sizes"] + " " + snakemake.output["bw"]

subprocess.call(bamCompare.split())
subprocess.call(wiggleTools.split())
subprocess.call(bdg2bw.split())
