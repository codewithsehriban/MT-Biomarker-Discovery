###load packages
library(Rsubread)
library(BiocParallel)
library(limma)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(dplyr)
library(data.table)

#featurecounts
# Set working directory to where your BAM files are located
setwd("/archive/d1/input/")
# Obtain BAM file path
bams <- list.files(pattern = ".bam", full.names = TRUE)

counts <- featureCounts(files = bams,
                          annot.ext = "/archive/mouse_gtf", 
                          isGTFAnnotationFile=TRUE,
                          GTF.featureType= c("exon"),
                          GTF.attrType="gene_id",
                          useMetaFeatures=TRUE,
                          allowMultiOverlap=F,
                          largestOverlap=F,
                          countMultiMappingReads=TRUE,
                          isPairedEnd=TRUE,
                          nthreads=10
)

counts_df <- as.data.frame(counts$counts)
write.table(counts$counts, "/archive/d1/counts.tsv", quote = F, col.names = T, sep = "\t")

#calculate fpkm values
z <- DGEList(counts=counts$counts, genes=counts$annotation[,c("GeneID","Length")])
z<- calcNormFactors(z)
RPKM <-rpkm(z)
write.table(RPKM, "/archive/d1/FPKM.tsv", quote = F, col.names = T, sep = "\t")
done
