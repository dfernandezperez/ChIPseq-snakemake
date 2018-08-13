#!/hpcnfs/software/r/3.5.0/bin/Rscript
# Usage --> input.bed output.txt distanceBeforeTSS distanceAfterTSS


####################                                                                                                              
# DEFINE FUNCTIONS #                                                                                              
####################

# Function to annotate the peaks using ChIPseeker
Peak_Annot <- function(infile, tssRegion = c(-2500, 2500)) {
  
  # Load packages and set parameters
  require(ChIPseeker)
  require(TxDb.Mmusculus.UCSC.mm9.knownGene)
  txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
  
  
  # read file and annotate peaks
  df <- readPeakFile(infile)
  annot <- annotatePeak(df, TxDb=txdb, annoDb="org.Mm.eg.db", tssRegion = tssRegion)
  
  # The program classifies promotres in 1kb, 2kb... this line removes that annotation and leaves just "Promoters"
  annot@anno$annotation <- sub(" \\(.*\\)", "", annot@anno$annotation)
  
  # The program changes the type of object, so go back to GRanges (useful for downstream analysis)
  final <- as.GRanges(annot)
  return(final)
}

args = commandArgs(trailingOnly=TRUE)

# Read inputs
file <- args[1]
outfile <- args[2]
before <- as.numeric(args[3])
after <- as.numeric(args[4])

annot <- Peak_Annot(file, tssRegion = c(-before, after))
write.table(as.data.frame(annot), file = outfile, sep = "\t", quote = F, row.names = F)

