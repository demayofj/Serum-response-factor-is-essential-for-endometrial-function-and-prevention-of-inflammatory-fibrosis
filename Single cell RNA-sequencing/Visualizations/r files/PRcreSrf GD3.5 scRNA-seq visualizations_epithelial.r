# Downstream analysis using object created by Sara Grimm
# Match the R version and important package versions used in the original analysis
# R432

library(Seurat);       # v5.0.1  
library(ggplot2);      # v3.4.4
library(dplyr);
library(viridis);
library(dittoSeq);
library(scales);
library(scCustomize)
library(tidyr)

# Record active package versions
packageVersion("Seurat");
packageVersion("ggplot2");        
packageVersion("dplyr");
packageVersion("viridis");
packageVersion("dittoSeq");
packageVersion("scales");
packageVersion("scCustomize");
packageVersion("tidyr");

# Read in the processed RDS file from Sara Grimm, which holds a seurat object with fibroblasts subsetted, integrated, and annotated (from 4 samples).
# We are using the broad clusters for these cell types
Adult_Srfflox_SrfKO_ep <- readRDS("/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-Epithelial.03feb2025.rds")
Adult_Srfflox_SrfKO_ep
head(Adult_Srfflox_SrfKO_ep)

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="broad_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

options(saved) # restore old settings

#Change colors

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE, 
        cols=c("0"="gold","2"="cornflowerblue","3"="darkorange1","4"="darkolivegreen2","1"="mediumorchid","5"="seagreen")
       )+NoLegend();

options(saved) # restore old settings

png(file = "Epi_PlotBySubcluster.png",
    width = 16, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE, 
        cols=c("0"="gold","2"="cornflowerblue","3"="darkorange1","4"="darkolivegreen2","1"="mediumorchid","5"="seagreen")
       )+NoLegend();

dev.off()

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE, 
        cols=c("0"="gold","2"="cornflowerblue","3"="darkorange1","4"="darkolivegreen2","1"="mediumorchid","5"="seagreen"),
        split.by="genotype"
        )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=FALSE, 
        cols=c("0"="gold","2"="cornflowerblue","3"="darkorange1","4"="darkolivegreen2","1"="mediumorchid","5"="seagreen"),
        split.by="genotype"
        )

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="EpiSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, 
        cols=c("GlandularEpithelial"="gold","ProliferatingGlandularEpithelial"="darkorange1","LuminalEpithelial"="mediumorchid","LowQualEpithelial"="seagreen"),
        )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="EpiSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, 
        cols=c("GlandularEpithelial"="gold","ProliferatingGlandularEpithelial"="darkorange1","LuminalEpithelial"="mediumorchid","LowQualEpithelial"="seagreen"),
        )

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="EpiSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, 
        cols=c("GlandularEpithelial"="gold","ProliferatingGlandularEpithelial"="darkorange1","LuminalEpithelial"="mediumorchid","LowQualEpithelial"="seagreen"),
        )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="genotype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE
       )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options


DimPlot(Adult_Srfflox_SrfKO_ep, group.by="EpiSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype",
        cols=c("GlandularEpithelial"="gold","ProliferatingGlandularEpithelial"="darkorange1","LuminalEpithelial"="mediumorchid","LowQualEpithelial"="seagreen"),
        )

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="EpiSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype",
        cols=c("GlandularEpithelial"="gold","ProliferatingGlandularEpithelial"="darkorange1","LuminalEpithelial"="mediumorchid","LowQualEpithelial"="seagreen"),
        )+NoLegend();

options(saved) # restore old settings

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$seurat_clusters

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="EpiSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype",
        cols=c("GlandularEpithelial"="gold","ProliferatingGlandularEpithelial"="darkorange1","LuminalEpithelial"="mediumorchid","LowQualEpithelial"="seagreen"),
        )+NoLegend();


options(saved) # restore old settings

# Take colors for use elsewhere
# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_ep_tmp <- Adult_Srfflox_SrfKO_ep
Idents(Adult_Srfflox_SrfKO_ep_tmp) <- Adult_Srfflox_SrfKO_ep_tmp$EpiSubtype
# Create vector with levels of object@ident

identities <- levels(Adult_Srfflox_SrfKO_ep_tmp)
# Create vector of default ggplot2 colors that can be used when creating plots
my_color_palette <- hue_pal()(length(identities))
my_color_palette

table(Adult_Srfflox_SrfKO_ep$EpiSubtype,Adult_Srfflox_SrfKO_ep$Sample);

table(Adult_Srfflox_SrfKO_ep$EpiSubtype,Adult_Srfflox_SrfKO_ep$genotype);

saved <- options(repr.plot.width = 8.5, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="Phase", pt.size=0.5, reduction='umap', shuffle=TRUE)

options(saved) # restore old settings

saved <- options(repr.plot.width = 15, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_ep, group.by="Phase", split.by="genotype", pt.size=0.5, reduction='umap', shuffle=TRUE)

options(saved) # restore old settings

table(Adult_Srfflox_SrfKO_ep$Phase,Adult_Srfflox_SrfKO_ep$genotype);

# Define marker genes

genes <- c("Epcam","Cdh1","Pgr","Esr1", #pan epithelial
           "Foxa2","Cxcl15","Sox9","Prss28","Prss29","Spink1","Ttr", #glandular
           "Calb1","Tacstd2","Lpar3", "Wnt7a", #luminal
           "Top2a","Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cdc20","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67","Tpx2", # high in all of G2M
           "Hells","Ung","Cdc6","Dtl", # S phase
           "Pclaf", "Ccna2" #Other proliferating (S/G2)
           );

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$seurat_clusters

saved <- options(repr.plot.width = 12, repr.plot.height = 3.5) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_ep, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

#Need to change ident level order to group clusters of the same cell type together
#Check ident order
Idents(Adult_Srfflox_SrfKO_ep)

#Copy to a temp file and change order of levels
Adult_Srfflox_SrfKO_ep_tmp <- Adult_Srfflox_SrfKO_ep

Idents(Adult_Srfflox_SrfKO_ep_tmp) <- Adult_Srfflox_SrfKO_ep_tmp$seurat_clusters

desired_order <- c("4", "1", "2", "0", "3", "5")
Idents(Adult_Srfflox_SrfKO_ep_tmp) <- factor(Idents(Adult_Srfflox_SrfKO_ep_tmp), levels = desired_order)
Idents(Adult_Srfflox_SrfKO_ep_tmp)

genes <- c("Epcam","Cdh1","Pgr","Esr1", #pan epithelial
           "Foxa2","Cxcl15","Sox9","Prss28","Prss29","Spink1","Ttr", #glandular
           "Calb1","Tacstd2","Lpar3", "Wnt7a", #luminal
           "Top2a","Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cdc20","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67","Tpx2", # high in all of G2M
           "Hells","Ung","Cdc6","Dtl", # S phase
           "Pclaf", "Ccna2" #Other proliferating (S/G2)
           );

saved <- options(repr.plot.width = 12, repr.plot.height = 3.5) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_ep_tmp, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

#Reduce list
genes <- c("Epcam","Cdh1","Pgr","Esr1", #pan epithelial
           "Foxa2","Cxcl15", #glandular
           "Calb1","Tacstd2","Lpar3", "Wnt7a", #luminal
           "Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cdc20","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67","Tpx2", # high in all of G2M
           "Hells","Ung","Cdc6","Dtl", # S phase
           "Pclaf", "Ccna2" #Other proliferating (S/G2)
           );

saved <- options(repr.plot.width = 12, repr.plot.height = 3.5) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_ep_tmp, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

saved <- options(repr.plot.width = 36, repr.plot.height = 24) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_ep_tmp, features = genes, 
            pt.size=0.5,
            label=FALSE,
            ncol=6
            )

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 72) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_ep_tmp, features = genes, 
            pt.size=0.5,
            label=FALSE,
            split.by="genotype",
            ncol=6
            )

options(saved) # restore old settings

#Reduce list more
genes <- c("Epcam","Cdh1","Pgr","Esr1", #pan epithelial
           "Foxa2","Cxcl15", "Sox9", #glandular
           "Calb1","Lpar3", "Tacstd2", #luminal
           "Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67", # high in all of G2M
           "Hells","Dtl" # S phase
           );

saved <- options(repr.plot.width = 7, repr.plot.height = 3.5) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_ep_tmp, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

# plot by subtype instead of cluster

genes <- c("Epcam","Cdh1","Pgr","Esr1", #pan epithelial
           "Foxa2","Cxcl15", "Sox9", #glandular
           "Calb1","Lpar3", "Tacstd2", #luminal
           "Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67", # high in all of G2M
           "Hells","Dtl" # S phase
           );

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 10, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_ep, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

# Split by genotype

genes <- c("Epcam","Cdh1","Pgr","Esr1", #pan epithelial
           "Foxa2","Cxcl15", "Sox9", #glandular
           "Calb1","Lpar3", "Tacstd2", #luminal
           "Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67", # high in all of G2M
           "Hells","Dtl" # S phase
           );

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 10, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_ep, split.by="genotype", cols=c("#6b00ff","#6b00ff"), features=genes) + RotatedAxis() + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

saved <- options(repr.plot.width = 36, repr.plot.height = 18) # Make the plots bigger from here on out, save old options

genes <- c("Epcam","Cdh1","Pgr","Esr1", #pan epithelial
           "Foxa2","Cxcl15", "Sox9", #glandular
           "Calb1","Lpar3", "Tacstd2", #luminal
           "Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67", # high in all of G2M
           "Hells","Dtl" # S phase
           );

FeaturePlot(Adult_Srfflox_SrfKO_ep, features = genes, 
            pt.size=0.5,
            label=FALSE,
            ncol=6
            )

options(saved) # restore old settings

# Use DittoSeq to plot fraction of each cluster by genotype


saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_ep, 
             var = "EpiSubtype", 
             group.by = "genotype",
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings             


# Use DittoSeq to plot fraction of each cluster by genotype (matching UMAP colors)


saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_ep, 
             var = "EpiSubtype", 
             group.by = "genotype",
             color.panel = c("gold",
                            "seagreen",
                            "mediumorchid",
                            "darkorange1"),
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings   

# Use DittoSeq to plot fraction of each cluster by genotype (matching UMAP colors) - turn


saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_ep, 
             var = "EpiSubtype", 
             group.by = "genotype",
             color.panel = c("gold",
                            "seagreen",
                            "mediumorchid",
                            "darkorange1"),
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20, angle = 90),
                    axis.text.y = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    title = element_text(size = 20, hjust = 0.5), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    aspect.ratio = 3/1)

# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Sprr2b","Sprr2d","Cxcl5","Lcn2","C3","Cxcl2","Ltf","Mmp7",
"Esr1","Egfr","Wnt4","Wnt7b","Hif1a","Jun","Muc1","Ezh2","Greb1","Clca3b",
"Ifnar1","Ifnar2","Ifngr1","Nfkb2","Tlr3","Tlr4","Irak4","Tnf","Rel","Relb","Hla-a","Ifi35","Irf1","Isg15","Isg20","Oas1a","Oas2","Oas3","Stat1","Stat2","Myd88","Nfkbia","Nfkbie","Tgfbr1","Zbp1","Nod1","S100a1","S100a9","S100a11",
"Il1r1","Il1rl2","Il33","Il6r","Il6st","Ccl27","Ccl28","Cxcl1","Cxcl14","Cxcl15","Cxcl16","Cxcl17","Cd44","Cd63","Mif",
"Pgr","Plcl1","Areg","Foxo1","Fkbp5",
"Ccnd1","Ccnd3","Cdc25a","Cdc34","Cdk4","Cdk6","Cdkn2b","E2f3","E2f5","Foxo1","Gnl3","Hdac2","Myc","Nrg1","Pa2g4","Smad3","Tfdp1","Tp53",
"Actb","Actg1","Ilk","Vim","Cnn2","Vit",
"Mmp14","Mmp23","Timp1","Timp3",
"Egr1",
"Casp3"
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
       )
options(saved) # restore old settings

# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 64) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Sprr2b","Sprr2d","Cxcl5","Lcn2","C3","Cxcl2","Ltf","Mmp7",
"Esr1","Egfr","Wnt4","Wnt7b","Hif1a","Jun","Muc1","Ezh2","Greb1","Clca3b",
"Ifnar1","Ifnar2","Ifngr1","Nfkb2","Tlr3","Tlr4","Irak4","Tnf","Rel","Relb","Ifi35","Irf1","Isg15","Isg20","Oas1a","Oas2","Oas3","Stat1","Stat2","Myd88","Nfkbia","Nfkbie","Tgfbr1","Zbp1","Nod1","S100a1","S100a9","S100a11",
"Il1r1","Il1rl2","Il33","Il6ra","Il6st","Ccl27a","Ccl28","Cxcl1","Cxcl14","Cxcl15","Cxcl16","Cxcl17","Cd44","Cd63","Mif",
"Pgr","Plcl1","Areg","Foxo1","Fkbp5",
"Ccnd1","Ccnd3","Cdc25a","Cdc34","Cdk4","Cdk6","Cdkn2b","E2f3","E2f5","Foxo1","Gnl3","Hdac2","Myc","Nrg1","Pa2g4","Smad3","Tfdp1","Trp53",
"Actb","Actg1","Ilk","Vim","Cnn2","Vit",
"Mmp14","Mmp23","Timp1","Timp3",
"Egr1",
"Casp3"
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4
       )

options(saved) # restore old settings

# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 56) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Sprr2b","Cxcl5","Lcn2","C3","Cxcl2","Ltf","Mmp7",
"Esr1","Egfr","Wnt4","Jun","Muc1","Ezh2","Greb1",
"Ifnar1","Ifnar2","Ifngr1","Nfkb2","Tlr3","Tlr4","Irak4","Tnf","Relb","Ifi35","Irf1","Isg15","Isg20","Oas1a","Oas2","Oas3","Stat1","Stat2","Myd88","Nfkbia","Tgfbr1","Zbp1","S100a1","S100a11",
"Il1r1","Il6st","Ccl28","Cxcl14","Cxcl15","Cxcl16","Cxcl17","Cd44","Cd63","Mif",
"Pgr","Plcl1","Areg","Foxo1","Fkbp5",
"Ccnd1","Ccnd3","Cdc25a","Cdc34","Cdk4","Cdk6","E2f3","E2f5","Foxo1","Gnl3","Hdac2","Nrg1","Pa2g4","Smad3","Tfdp1","Trp53",
"Actb","Actg1","Ilk",
"Mmp14","Mmp23","Timp3",
"Egr1",
"Casp3"
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Remove low expressors, remove some labels, add dots
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 56) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Sprr2b","Cxcl5","Lcn2","C3","Cxcl2","Ltf","Mmp7",
"Esr1","Egfr","Wnt4","Jun","Muc1","Ezh2","Greb1",
"Ifnar1","Ifnar2","Ifngr1","Nfkb2","Tlr3","Tlr4","Irak4","Tnf","Relb","Ifi35","Irf1","Isg15","Isg20","Oas1a","Oas2","Oas3","Stat1","Stat2","Myd88","Nfkbia","Tgfbr1","Zbp1","S100a1","S100a11",
"Il1r1","Il6st","Ccl28","Cxcl14","Cxcl15","Cxcl16","Cxcl17","Cd44","Cd63","Mif",
"Pgr","Plcl1","Areg","Foxo1","Fkbp5",
"Ccnd1","Ccnd3","Cdc25a","Cdc34","Cdk4","Cdk6","E2f3","E2f5","Foxo1","Gnl3","Hdac2","Nrg1","Pa2g4","Smad3","Tfdp1","Trp53",
"Actb","Actg1","Ilk",
"Mmp14","Mmp23","Timp3",
"Egr1",
"Casp3"
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0.001,
        ncol = 4
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Reduce list
# Remove low expressors, remove some labels, add dots
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Sprr2b","Lcn2","C3","Ltf","S100a1", #epithelial innate immunity
"Cxcl2","Cxcl5","Cxcl15","Cxcl17","Mif", #cytokines
"Mmp7","Mmp14","Timp3", #mmps
"Wnt4","Muc1","Greb1", #E2 signaling
"Nfkb2","Relb","Tlr4","Irak4", #Nfkb signaling
"Ifi35","Irf1","Isg15","Oas1a","Oas2","Stat1","Stat2", #interferon signaling
"Zbp1","Casp3", #cell death process
"Ccnd1","Ccnd3","Cdc25a","Cdc34","Cdk4","Cdk6", #cell cycle
"Actb","Actg1", #cytoskeleton
"Egr1","Egfr", #growth response
"Plcl1","Areg","Fkbp5" #progesterone signaling
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0.001,
        ncol = 4
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Reduce list more, change layout
# Remove low expressors, remove some labels, add dots
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Lcn2","C3","Ltf", "Zbp1", #epithelial innate immunity
"Cxcl5","Cxcl15","Cxcl17","Mif", #cytokines
"Mmp7","Mmp14","Timp3", #mmps
"Wnt4","Muc1","Greb1", #E2 signaling
"Nfkb2","Relb","Tlr4","Irak4", #Nfkb signaling
"Irf1","Isg15","Oas1a","Stat1", #interferon signaling
"Casp3" #apoptosis
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0.001,
        ncol = 6
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Reduce list more, change layout, remove dots
# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Lcn2","C3","Ltf", "Zbp1", "S100a1", "Mif", #epithelial innate immunity
"Cxcl5","Cxcl15","Cxcl17", #chemokines
"Irf1","Isg15","Stat1", #interferon signaling
"Mmp7","Mmp14","Timp3", #mmps
"Greb1","Wnt4","Muc1" #E2 signaling
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Alternate
# Reduce list more, change layout, remove dots
# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Lcn2","C3","Ltf", "Zbp1", "Mif", #epithelial innate immunity
"Egr1", #growth response
"Cxcl5","Cxcl15","Cxcl17", #chemokines
"Irf1","Isg15","Stat1", #interferon signaling
"Mmp7","Mmp14","Timp3", #mmps
"Greb1","Wnt4","Muc1" #E2 signaling
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# More from human
# Reduce list more, change layout, remove dots
# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 2) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Dcn","Gata2","Muc4", "Csf1", "Fbln1"
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# more as a follow up
# Reduce list more, change layout, remove dots
# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 2) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Rara","Rbp1", "Foxo1", "Pgr"
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Alternate update
# Reduce list more, change layout, remove dots
# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Lcn2","C3","Ltf",  #epithelial innate immunity
"Mmp7","Timp3","Fbln1", #ECM
"Cxcl15","Cxcl17","Mif", #cytokines
"Isg15","Stat1", "Oas3", #interferon signaling
"Greb1","Muc4","Wnt4", #E2 signaling
"Gata2", "Fkbp5", #Progesterone signaling
"Egr1" #growth response
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Another alternate
# Reduce list more, change layout, remove dots
# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Lcn2","C3","Ltf",  #epithelial innate immunity
"Mmp7","Timp3","Fbln1", #ECM
"Egr1", #growth response
"Cxcl17","Mif", #cytokines
"Isg15","Stat1","Oas3", #interferon signaling
"Actb", "Tpm1", #cytoskeleton
"Esr1","Greb1","Muc4", #E2 signaling
"Pgr", "Gata2", "Fkbp5", "Areg" #Progesterone signaling
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings        


# Another alternate
# Reduce list more, change layout, remove dots
# Remove low expressors, remove some labels
# Violin Plot luminal DEGs

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Esr1","Greb1","Muc4", #E2 signaling
"Lcn2","Ltf","C3","Mif",  #innate immunity
"Cxcl17","Mmp7","Timp3","Fbln1", #ECM and fibrosis
"Isg15","Stat1","Oas3", #interferon signaling
"Pgr", "Gata2", "Fkbp5", "Areg", #Progesterone signaling
"Actb", "Tpm1", #cytoskeleton
"Egr1" #growth response
), 
        idents = "LuminalEpithelial", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings        


# Violin Plot glandular DEGs and markers

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 56) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Egr1","Egr2","Egr3","Egr4","Egfr","Fos","Jun","Fosb","Stat3", "Hbegf",#growth response
"Foxa2","Prss29","Spink3","Cxcl15","Lif","Ttr","Spink1","Prss23", #gland markers (Foxa2,Prss29,Spink3,Cxcl15 not DEGs)
"Areg","Fkbp5","Foxo1","Pgr", #P4 signaling
"Vim","Itga5","Actb","Lpar1","Rhoj","Acta2","Actg1","Myl9", "Cnn2", "Tpm1", #cytoskeleton
"Muc4","Greb1", "Esr1", #estrogen signaling (Esr1 is not DEG)
"Cenpf","Mki67","E2f2","Cdk1","Hells", #proliferation
"Mcm3","Mcm4","Mcm6","Mcm7", #cell cycle
"Cxcl17","Tnfsf10","Il10rb","Mif", #cytokines
"Mmp7","Fbln1","Timp3","Timp4", #ECM
"Ifit3","Oas2","Isg15","Ifit1","Ifi44","Ifi47","Irf7","Ifi27","Oasl1","Stat1","Ifitm1","Stat2","Oas1a","Oas1g","Tlr3","Oasl2","Irf9","Irf1", #interferon
"Lcn2","Zbp1","C3","S100a9","Ltf","Sprr2f","C2", #innate inflammation
"Fbln7","Dio2","Cxcl14","Cebpb", #decidualization/fibroblast
"Elf5", #E2 in breast cancer
"Calb1","Cdh1","Krt19","Epcam", #epithelial
"Rara", #retinoic acid signaling
"Trp53", #p53
"Ptgs1" #Cox1
), 
        idents = "GlandularEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Remove low expressors/uninteresting
# Violin Plot glandular DEGs and markers

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Egr1","Egr3","Egr4","Egfr","Fos","Jun","Fosb","Stat3", "Hbegf", #growth response
"Foxa2","Cxcl15","Lif","Ttr","Spink1", #gland markers (Foxa2,Prss29,Spink3,Cxcl15 not DEGs)
"Areg","Fkbp5","Foxo1","Pgr", #P4 signaling
"Itga5","Actb","Lpar1","Rhoj","Actb","Actg1", "Tpm1", #cytoskeleton
"Muc4","Greb1", "Esr1", #estrogen signaling (Esr1 is not DEG)
"Cxcl17","Mif", #cytokines
"Mmp7","Fbln1","Timp3", #ECM
"Isg15","Ifit1","Ifi44","Ifi47","Irf7","Ifi27","Oasl1","Stat1","Ifitm1","Stat2","Oas1a","Tlr3","Oasl2","Irf9","Irf1", #interferon
"Lcn2","Zbp1","C3","Ltf","Sprr2f","C2" #innate inflammation
), 
        idents = "GlandularEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Reduce list
# Violin Plot glandular DEGs and markers

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 16) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Egr1","Egfr","Hbegf", #growth response
"Foxa2","Cxcl15","Ttr","Spink1","Lif", #gland markers (Foxa2,Prss29,Spink3,Cxcl15 not DEGs)
"Pgr","Areg","Fkbp5", #P4 signaling
"Esr1","Muc4","Greb1", #estrogen signaling (Esr1 is not DEG)
"Cxcl17","Mif", #cytokines
"Mmp7","Fbln1","Timp3", #ECM
"Actb","Actg1", "Tpm1", #cytoskeleton
"Isg15","Irf7","Ifi27","Oasl1","Stat1","Ifitm1", "Oas1a","Oasl2","Irf9","Irf1", #interferon
"Lcn2","Zbp1","C3","Ltf","Sprr2f" #innate inflammation
), 
        idents = "GlandularEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Reduce list more, change layout
# Violin Plot glandular DEGs and markers

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Egr1","Egfr","Hbegf", #growth response
"Foxa2","Cxcl15","Ttr","Spink1","Lif", #gland markers (Foxa2,Prss29,Spink3,Cxcl15 not DEGs)
"Pgr","Areg","Fkbp5", #P4 signaling
"Esr1","Muc4","Greb1", #estrogen signaling (Esr1 is not DEG)
"Cxcl17","Mif", #cytokines
"Mmp7","Fbln1","Timp3", #ECM
"Isg15","Stat1","Irf7","Oas1a","Oasl2", #interferon
"Lcn2","C3","Ltf","Sprr2f", #innate inflammation
"Actb","Actg1", "Tpm1" #cytoskeleton
), 
        idents = "GlandularEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 7
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Reduce list more, change layout
# Violin Plot glandular DEGs and markers

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Lcn2","C3","Ltf","Sprr2f", #innate inflammation
"Mmp7","Fbln1","Timp3", #ECM
"Cxcl17","Mif", #cytokines
"Isg15","Stat1","Oasl2", #interferon
"Egr1","Hbegf", #growth response
"Esr1","Muc4","Greb1", #estrogen signaling (Esr1 is not DEG)
"Pgr","Areg","Fkbp5", #P4 signaling
"Foxa2","Cxcl15","Ttr","Spink1","Lif", #gland markers (Foxa2,Prss29,Spink3,Cxcl15 not DEGs)
"Actb","Actg1", "Tpm1" #cytoskeleton
), 
        idents = "GlandularEpithelial", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 7
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings

# Reduce list more, change layout, update colors
# Violin Plot glandular DEGs and markers

Idents(Adult_Srfflox_SrfKO_ep) <- Adult_Srfflox_SrfKO_ep$EpiSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_ep, features = c(
"Esr1","Greb1","Muc4", #estrogen signaling (Esr1 is not DEG)
"Lcn2","Ltf","Sprr2f","Mif", #innate inflammation
"Cxcl17","Mmp7","Timp3","Fbln1",#ECM
"Isg15","Stat1","Oasl2", #interferon
"Pgr","Areg","Fkbp5", #P4 signaling
"Actb","Tpm1", #cytoskeleton
"Egr1","Hbegf", #growth response
"Foxa2","Cxcl15","Ttr","Spink1","Lif" #gland markers (Foxa2,Prss29,Spink3,Cxcl15 not DEGs)
), 
        idents = "GlandularEpithelial", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7
       ) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

options(saved) # restore old settings
