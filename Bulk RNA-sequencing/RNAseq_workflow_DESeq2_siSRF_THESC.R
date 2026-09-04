# Ryan Marquardt
# Modified from Updated_RNAseq_DESeq2_workflow_2020 by Jake Reske and EdgeR scrpt from Ty Wang, 2022; converted to DESeq2 with help from Peter Lais

# RNA-seq experimental design n=3 replicates/treatment siNT vs siSRF, Veh only.
# Cells: THESC.
# Goal: use DESeq2 to perform differential gene expression analysis beginning with a raw counts table/matrix of m genes by n samples.
###############################################
###############################################

# load dependencies (install first if necessary-google it)
library("ggplot2")
library("DESeq2")
library("dplyr")

###############################################
# set working directory
setwd("//wine/DeMayo_Francesco Group/Xdata/Ryan Marquardt/1_Projects/1_SRF/My Bioinformatics/siSRF THESC RNA-seq/DESeq2/250315 test with veh only")

# import raw counts data
raw.counts <- read.csv("feature_counts_import_250315.csv")
# convert to numeric matrix
ids <- raw.counts$Gene
raw.counts <- raw.counts[, -1] # remove gene ID column
raw.counts <- as.matrix(raw.counts)
rownames(raw.counts) <- ids
rm(ids)

###############################################

# prepare samples table containing identifiers
samples <- data.frame(c(colnames(raw.counts)))
# rename samples column
colnames(samples)[1] <- "sample"
# add sample groups
samples$Treatment <- factor(c("NT","NT","NT","SRF","SRF","SRF"))

# add row names
rownames(samples) <- samples$sample

# relevel to set Control as reference, or another level if desired
samples$Treatment <- relevel(samples$Treatment, "NT")

###############################
# DESeq2 setup for combined PCA
dds <- DESeqDataSetFromMatrix(countData = raw.counts,
                              colData = samples,
                              design = ~ Treatment)
#output was Large DESeqDataSet (59251 elements, 9.8 MB)

# Pre-filter low count genes (minimum set to average of 1 count per sample)
keep <- rowSums(counts(dds)) >= ncol(counts(dds))
dds <- dds[keep, ]
rm(keep)
# 27,185 expressed gene remain

##################
# calculate normalized counts (median of ratios method)
dds <- estimateSizeFactors(dds)
normalized.counts <- counts(dds, normalized=TRUE)

# write normalized counts output to csv
write.csv(normalized.counts, file="normalized_counts.csv")

#################
# regularized-logarithm transformation (rlog) for downstream PCA/MDS etc.
# note: blind=FALSE elicits slightly-supervised transformation
# "differences between [experimental design variables] will not contribute to expected variance-mean trend of data"
# blind=TRUE elicits fully unsupervised transformation
rld <- rlog(dds, blind = FALSE)
rlog.counts <- assay(rld)

# write rlog counts output to csv
write.csv(rlog.counts, file="rlog_counts.csv")

##########################################

# plot PCA based on rlog-transformed counts

#PCA per sample:
pcaData.rlog <- plotPCA(rld, intgroup=c("sample"), returnData=TRUE)

percentVar.rlog <- round(100 *attr(pcaData.rlog, "percentVar"))

ggplot(pcaData.rlog, aes(x = PC1, y = PC2, color = sample)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar.rlog[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar.rlog[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

#PCA per treatment:
pcaData.rlog <- plotPCA(rld, intgroup=c("Treatment"), returnData=TRUE)

percentVar.rlog <- round(100 *attr(pcaData.rlog, "percentVar"))

ggplot(pcaData.rlog, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar.rlog[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar.rlog[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())



#####################	DEG list ####################################


# Perform the DESeq2 normalization and get the results.
dds <- DESeq(dds)
res <- results(dds)

# Convert the results to a data.frame and save. You can later filter this in Excel.
res_df <- as.data.frame(res)
write.table(res_df, sep = "\t", file = "DEGs-NTveh-SRFveh_250315.txt");
### Note in output file it will be necessary to shift the column headers over by one cell

# As a preview, print the number of genes with adjusted P-value < 0.05.
res_df_no_na <- res_df[!is.na(res_df$padj),]
sum(res_df_no_na$padj < 0.05)
# 4062

