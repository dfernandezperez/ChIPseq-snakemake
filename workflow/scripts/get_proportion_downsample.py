import pysam
import argparse

##-----------------------------------------------------------------
## ARGUMENT PARSING
##-----------------------------------------------------------------
parser = argparse.ArgumentParser(description='Calculate the target number of reads to downsample a list of bam files.')
parser.add_argument('-i', '--input', help='List of input bam files', action='append', required = True, nargs='+')
parser.add_argument('-o', '--output', help='outpur file',  required = True)

options = parser.parse_args()

output    = options.output
inputs    = options.input
inputs    = inputs[0]

##-----------------------------------------------------------------
## Code
##-----------------------------------------------------------------
nreads = []

for i in inputs:
    nreads.append(pysam.AlignmentFile(i, "rb").mapped)

min_reads = min(nreads)

with open(output, 'w') as f:
    f.write(str(min_reads))