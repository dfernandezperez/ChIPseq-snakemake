import subprocess
import pysam

def touch_file(bam):
	command = "touch " + bam
	subprocess.call(command.split())

# Avoid the warning from pysam that tells that the bam file is more recent than the index file (due to snakemake behaviour)
touch_file(snakemake.input["case"] + ".bai")

####################
## Read bam files ##
####################
reference = pysam.AlignmentFile(snakemake.input["reference"], "rb")
case = pysam.AlignmentFile(snakemake.input["case"], "rb")

################################################################
## Calculate normalization factors: (1/mapped reads)*1million ##
################################################################
reference_norm = str( (1.0/float(reference.mapped))*1000000 )
case_norm = str( (1.0/float(case.mapped))*1000000 )
print("Scaling factor for " + snakemake.input["case"] + " is " + case_norm)

#############################
## Bash commands to launch ##
#############################
bamCompare = "bamCompare -b1 " + snakemake.input["case"] + " -b2 " + snakemake.input["reference"] + " -o " + snakemake.output["bw"] + " -bs 50 --ratio subtract --scaleFactors " + case_norm + ":" + reference_norm + " " + snakemake.threads + " --extendReads 200 &> " + snakemake.log
subprocess.call(bamCompare.split())