# Serum-response-factor-is-essential-for-endometrial-function-and-prevention-of-inflammatory-fibrosis
Code associated with the following publication:
Marquardt RM, Grimm SA, Wu SP, Lais PF, Li SY, Xu X, Smithberger E, Cunefare D, Ganta C, Olson D, Jeong EM, Jeong JW, Lessey BA, Lydon JP, DeMayo FJ. Serum response factor is essential for endometrial function and prevention of inflammatory fibrosis. Proc Natl Acad Sci U S A. 2025 Nov 4;122(44):e2510060122. doi: 10.1073/pnas.2510060122. Epub 2025 Oct 28. PMID: 41150713; PMCID: PMC12595411.

## Single cell RNA-sequencing
The following command-line tools and R packages were used in preparing the results:
**10x Genomics Cell Ranger, scDblFinder, SoupX, Seurat, dittoSeq, CellChat, MuDataSeurat**
Version numbers for each tool is specified in the manuscript supporting material.

See the `Single cell RNA-Sequencing` folder for code used to process and analyze the datasets.

## Bulk RNA-sequencing

The following external tools were used in preparing the results:

* Command-line tools: **Trimmomatic, HISAT2, SAMTools, Subread’s featureCounts tool**
* R packages: **DESeq2, ComplexHeatmap, circlize**

Version numbers for each tool is specified in the manuscript supporting material. For details on use of the command-line tools (used to create count matrices), please see the [DeMayo Bulk Omics](https://github.com/Petronian/demayo-bulk-omics) repository.

See the `Bulk-RNA-Sequencing` folder for code used to create and analyze DEG lists for the bulk RNA-seq datasets.
