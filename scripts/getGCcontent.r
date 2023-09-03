#!/usr/bin/Rscript
# +
## Devon Ryan - dpryan79
# -
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

GTFfile = "/research/projects/chlamydomonas/MAexpression/data/genome_info/v6_genome_plus_anno/CC4532.v1_1.gene_exons.gtf"
FASTAfile = "/research/projects/chlamydomonas/MAexpression/data/genome_info/v6_genome_plus_anno/CC4532.w_organelles_MTplus.fa"

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="CC4532.v1_1", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementLengths=elementNROWS
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementLengths(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")

write.table(output, file="/research/projects/chlamydomonas/MAexpression/data/genome_info/v6_genome_plus_anno/GC_lengths.tsv", sep="\t")
