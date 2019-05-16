import pysam
import sys

def get_all_readsId(bamfile):
# Funtion to get a list (set) with all the read IDs from a bamfile. First it will index the bam.

	f = pysam.AlignmentFile(bamfile,'rb')
	all_readsid = set()
	for read in f.fetch():
		all_readsid.add(read.qname)

	f.close()
	return(all_readsid)

def get_common_readId(bamfile1, bamfile2):
# Get the common elements from 2 sets (wich contain the read ids from each bam)

	ids_1 = get_all_readsId(bam1)
	ids_2 = get_all_readsId(bam2)

	common = ids_1.intersection(ids_2)

	return(common)

def write_unique_reads(bamfile1, bamfile2, common_reads):
# This function will open the bam files and write to a new ones just those reads that aren't shared by both bam files

	bam1_open = pysam.AlignmentFile(bamfile1,'rb')
	bam2_open = pysam.AlignmentFile(bamfile2,'rb')

	# Open the new bam files where were going to write the output
	bam1_outfile = pysam.Samfile(bamfile1 + '.temporary', 'wb', template=bam1_open)
	bam2_outfile = pysam.Samfile(bamfile2 + '.temporary', 'wb', template=bam2_open)

	# Print the bam lines that are not shared by mm and dm to new bam files
	counter = 0
	for read in bam1_open.fetch():
	    if read.qname not in common_reads:
	        bam1_outfile.write(read)

	for read in bam2_open.fetch():
	    if read.qname not in common_reads:
	        bam2_outfile.write(read)
	    else:
	        counter +=1

	print "Removed %s reads" % (counter)



if __name__ == "__main__":
	bam1, bam2 = sys.argv[1], sys.argv[2]
	common = get_common_readId(bam1, bam2)
	write_unique_reads(bam1, bam2, common)
