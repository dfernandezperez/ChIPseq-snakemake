import subprocess
import pysam
import os
import argparse

######################
## ARGUMENT PARSING ##
######################
parser = argparse.ArgumentParser(description='Bam to bigwig for spike-in samples')
parser.add_argument('-c', '--case', help='case sample bam file', required=True)
parser.add_argument('-s', '--spike', help='spike-in bam file', required=True)
parser.add_argument('-i', '--refmm', help='reference input bam with spike mapped to reference', required=True)
parser.add_argument('-x', '--refdm', help='reference input bam with spike mapped to spike', required=True)
parser.add_argument('-b', '--bigwig', help='name for the output bigwig file', required=True)
parser.add_argument('-p', '--threads', help='Number of threads to use', required=True)
parser.add_argument('-o', '--otherParams', help='Extra parameters to deeptools', action='append', nargs=argparse.REMAINDER)

options = parser.parse_args()

params    = options.otherParams
params    = " ".join(str(e) for e in params[0])
case      = options.case
spike     = options.spike
ref_input = options.refmm
ref_spike = options.refdm
threads   = options.threads
bw        = options.bigwig


##############################
## Avoid warning from pysam ##
##############################
def touch_file(file):
	command = "touch " + file
	subprocess.call(command.split())

# Avoid the warning from pysam that tells that the bam file is more recent than the index file (due to snakemake behaviour)
touch_file(spike + ".bai")
touch_file(case + ".bai")

####################
## Read bam files ##
####################
dm        = pysam.AlignmentFile(spike, "rb")
c         = pysam.AlignmentFile(case, "rb")
ref_input = pysam.AlignmentFile(ref_input, "rb")
ref_spike = pysam.AlignmentFile(ref_spike, "rb")

################################################################
## Calculate normalization factors: (1/mapped reads)*1million ##
################################################################
# Following https://bio-protocol.org/e2981#biaoti24681
Nm    = dm.mapped
gamma = float(ref_spike.mapped)/float(ref_input.mapped)
alfa = str(gamma/Nm*1000000)
print alfa

#############################
## Bash commands to launch ##
#############################
bamCoverage = "bamCoverage -b " + case + " -o " + bw + " --scaleFactor " + alfa + " -p " + threads + " " + params

subprocess.call(bamCoverage.split())
