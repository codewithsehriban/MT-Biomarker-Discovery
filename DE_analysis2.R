library(Rsubread)
library(BiocParallel)
library(limma)
library(edgeR)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(dplyr)
library(data.table)

# Set working directory
setwd("/input/")

# Get file paths for BAM files in directories
bams_S <- list.files(path = ".", pattern = ".bam", full.names = TRUE)

# Run featureCounts for single-end data
counts_S <- featureCounts(files = bams_S,
                          annot.ext = "/archive/mouse_gtf", 
                          isGTFAnnotationFile = TRUE,
                          GTF.featureType = "exon",
                          GTF.attrType = "gene_id",
                          useMetaFeatures = TRUE,
                          allowMultiOverlap = FALSE,
                          largestOverlap = FALSE,
                          countMultiMappingReads = TRUE,
                          isPairedEnd = TRUE,
                          nthreads = 10)

# Extract counts and annotations for single-end data
counts_S_df_C <- as.data.frame(counts_S$counts)
counts_S_df_A <- as.data.frame(counts_S$annotation)
counts_S_df_C$GeneID <- row.names(counts_S_df_C)
counts_S_df_C_A <- merge(counts_S_df_A, counts_S_df_C, by = "GeneID")

# Import annotation file
gtf <- rtracklayer::import("/archive/mouse_gtf")
gtf_data_frame <- as.data.frame(gtf)

# Filter annotation data
gtf_filtered <- select(gtf_data_frame, c(10, 12, 17))
gtf_filtered_uniq <- gtf_filtered %>% distinct(gene_id, .keep_all = TRUE)

# Add gene IDs as a column
counts_S_df_C_A$gene_id <- row.names(counts_S_df_C_A)

# Merge counts and annotations
counts_df_annotated <- merge(gtf_filtered_uniq, counts_S_df_C_A, by = "gene_id")

# Write annotated counts to file
write.table(counts_df_annotated, "counts_df_annotated.tsv", quote = FALSE, col.names = TRUE, sep = "\t")

# Calculate FPKM values
z <- DGEList(counts = counts_S$counts, genes = counts_S_df_C_A[, c("GeneID", "Length")])
n <- calcNormFactors(z)
RPKM <- rpkm(n)
fp <- as.data.frame(RPKM)
fp$gene_id <- rownames(fp)
FPKM <- fp %>% select(gene_id, everything())
FPKM_annotated <- merge(FPKM, gtf_filtered_uniq, by = "gene_id")

# Write FPKM values to file
write.table(FPKM_annotated, "pkm_annotated.tsv", quote = FALSE, col.names = TRUE, sep = "\t")
