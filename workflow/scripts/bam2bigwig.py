import subprocess
import pysam
import argparse

######################
## ARGUMENT PARSING ##
######################
parser = argparse.ArgumentParser(description='Bam to bigwig with deeptools')
parser.add_argument('-c', '--case', help='case sample bam file', required=True)
parser.add_argument('-r', '--reference', help='reference file (input) bam file', required=True)
parser.add_argument('-b', '--bigwig', help='name for the output bigwig file', required=True)
parser.add_argument('-p', '--threads', help='Number of threads to use', required=True)
parser.add_argument('-o', '--otherParams', help='Extra parameters to deeptools', action='append', nargs=argparse.REMAINDER)

options = parser.parse_args()

params    = options.otherParams
params    = " ".join(str(e) for e in params[0])
case      = options.case
reference = options.reference
threads   = options.threads
bw        = options.bigwig

##############################
## Avoid warning from pysam ##
##############################
def touch_file(file):
	command = "touch " + file
	subprocess.call(command.split())

# Avoid the warning from pysam that tells that the bam file is more recent than the index file (due to snakemake behaviour)
touch_file(case + ".bai")

####################
## Read bam files ##
####################
rf = pysam.AlignmentFile(reference, "rb")
c  = pysam.AlignmentFile(case, "rb")

################################################################
## Calculate normalization factors: (1/mapped reads)*1million ##
################################################################
reference_norm = str( (1.0/float(rf.mapped))*1000000 )
case_norm      = str( (1.0/float(c.mapped))*1000000 )


#############################
## Bash commands to launch ##
#############################
bamCompare = "bamCompare -b1 " + case + " -b2 " + reference + " -o " + bw + " --operation subtract --scaleFactors " + case_norm + ":" + reference_norm + " -p " + threads + " " + params

subprocess.call(bamCompare.split())
