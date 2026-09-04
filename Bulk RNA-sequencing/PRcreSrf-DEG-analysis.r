library(DESeq2)
library(pheatmap)
library(data.table)
library(purrr)
library(tidyr)
library(dplyr)
library(Seurat)
library(scatterplot3d)
library(biomaRt)
library(ggplot2)
library(cowplot)

getwd()

# Load in the table and make sure it works/reads fine!
# There's a comment at the beginning of the counts file, so skip that first line.
filePath = "analysis/counts.txt"
cts = fread(file=filePath)
head(cts)

# Clean up the columns of the counts matrix.
cleanedColNames = unlist(map(strsplit(basename(colnames(cts)[7:ncol(cts)]), '_'), function (x) x[[1]]))
colnames(cts) = c(colnames(cts)[1:6], cleanedColNames)
head(cts)

# Make the DESeq counts matrix.
dcts = as.data.frame(cts[,7:ncol(cts)])
rownames(dcts) = unlist(cts[,'Geneid'])
head(dcts)

# Make a factor list between control and experimental groups.
# We can do this by just stripping the numbers out of the cleaned group names.
groups = as.factor(sub("\\d", "", cleanedColNames))
groups = relevel(groups, ref="Ctrl")
dgroups = data.frame(row.names=cleanedColNames, Group=groups)
dgroups
levels(groups)

# Ensure the data and metadata labels are within each other and in the same order.
all(rownames(groups) %in% colnames(dcts))
all(rownames(groups) == colnames(dcts))

# Normalize and center the data in advance of running PCA.
dctsProc = NormalizeData(dcts)
dctsProc = ScaleData(dctsProc)
dctsProc = as.matrix(dctsProc)
head(dctsProc, 8)

# Calculate principal components. Note that things have already been centered and scaled.
pca_results = prcomp(t(dctsProc))
pca_results$x

# Make a scatterplot with the 3D PCA data.
# To help pinpoint them in 3D space, I project them down to the XY plane as well.
colors_reference = c("#56B4E9", "#E69F00")
legend = c("Ctrl", "Srf")
color_indices = as.numeric(startsWith(rownames(pca_results$x[,1:3]), 'Srf')) + 1
colors = colors_reference[color_indices]

plt = scatterplot3d(pca_results$x[,1:3], color=colors, pch=16, type="h")

legend(plt$xyz.convert(210, 0, 50), legend=legend, col=colors_reference, pch = 16)

# Make the DESeqDataSet and analyze.
dds = DESeqDataSetFromMatrix(countData = dcts,
                             colData = dgroups,
                             design = ~ Group)
dds

dds = DESeq(dds)
res = results(dds)
res

# Carry out the above instructions.
testCounts = counts(dds)
testCounts = cbind(testCounts, res[['padj']], 2 ** res[['log2FoldChange']])
colnames(testCounts) = c(colnames(testCounts)[1:8], 'padj', 'FC')
testCountsdf = as.data.frame(testCounts)
testCountsdf = testCountsdf[!is.na(testCountsdf[['FC']]),]
head(testCountsdf)

temp = testCountsdf[['FC']]
temp[temp < 1] = -1 / temp[temp < 1]
testCountsdf[['FC']] = temp
head(testCountsdf)

testCountsThresholdeddf = subset(testCountsdf, padj < 0.05 & abs(FC) > 1.5)
testCountsThresholdeddf = testCountsThresholdeddf[order(testCountsThresholdeddf[['FC']], decreasing=TRUE),]

# Add in MGI symbol column. We used ensembl version 102 for this.
ensembl = useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl", version=102)

mapperTable = subset(getBM(c('ensembl_gene_id', 'mgi_symbol'), mart=ensembl), mgi_symbol!="")
mapper = mapperTable[['mgi_symbol']]
names(mapper) = mapperTable[['ensembl_gene_id']]

geneMapper = function(ensembl_id) {
    if (!(ensembl_id %in% names(mapper))) return(ensembl_id)
    else return(mapper[ensembl_id])
}

testCountsThresholdeddf[['MGI-Symbol']] = map(rownames(testCountsThresholdeddf), geneMapper)

head(testCountsThresholdeddf)
tail(testCountsThresholdeddf)
nrow(testCountsThresholdeddf)

# Output results to a CSV
fwrite(testCountsThresholdeddf, "./analysis/topDEGs-posNeg.csv", row.names=TRUE)
