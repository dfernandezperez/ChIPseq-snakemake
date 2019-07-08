import subprocess
import pysam
from argparse import ArgumentParser

######################
## ARGUMENT PARSING ##
######################
parser = ArgumentParser(description='Bam to bigwig with deeptools')
parser.add_argument('-c', '--case', help='case sample bam file', required=True)
parser.add_argument('-b', '--bigwig', help='name for the output bigwig file', required=True)
parser.add_argument('-p', '--threads', help='Number of threads to use', required=True)

options = parser.parse_args()

case      = options.case
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
c  = pysam.AlignmentFile(case, "rb")

################################################################
## Calculate normalization factors: (1/mapped reads)*1million ##
################################################################
case_norm      = str( (1.0/float(c.mapped))*1000000 )


#############################
## Bash commands to launch ##
#############################
bamCoverage = "bamCoverage -b " + case + " -o " + bw + " -e -p " + threads + " --scaleFactor " + case_norm

subprocess.call(bamCoverage.split())
