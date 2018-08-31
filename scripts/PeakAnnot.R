# Usage --> input.bed output.txt distanceBeforeTSS distanceAfterTSS
require(dplyr)


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


#################
## Read inputs ##
#################
input <- args[1]
before <- as.numeric(args[2])
after <- as.numeric(args[3])
out1 <- args[4]
out2 <- args[5]
out3 <- args[6]
out4 <- args[7]
out5 <- args[8]

####################
## Annot bed file ##
####################
annot <- Peak_Annot(input, tssRegion = c(-before, after)) %>% as.data.frame()

# Filter annotared bed file to obtain promoter target genes, bed files from peaks overlaping with promoters/distal regions and the coordinates of the promoter target genes.
distal.peaks <- annot %>% subset(annotation != "Promoter") %>% dplyr::select(c("seqnames", "start", "end", "V4", "V5")) 
promo.peaks <- annot %>% subset(annotation == "Promoter") %>% dplyr::select(c("seqnames", "start", "end", "V4", "V5")) 
promo.targets.bed <-  annot %>% subset(annotation == "Promoter") %>% 
  dplyr::select(c("seqnames", "geneStart", "geneEnd", "ENSEMBL", "SYMBOL", "geneStrand")) %>% 
  mutate(geneStrand = replace(geneStrand, geneStrand == 1, "+")) %>%
  mutate(geneStrand = replace(geneStrand, geneStrand == 2, "-")) %>% unique()
promo.targets <- annot %>% subset(annotation == "Promoter") %>% dplyr::select("SYMBOL") %>% unique()


##################
## Write output ##
##################
write.table(annot, file = out1, sep = "\t", quote = F, row.names = F)
write.table(promo.targets.bed, file = out2, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(promo.targets, file = out3, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(promo.peaks, file = out4, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(distal.peaks, file = out5, sep = "\t", quote = F, row.names = F, col.names = F)