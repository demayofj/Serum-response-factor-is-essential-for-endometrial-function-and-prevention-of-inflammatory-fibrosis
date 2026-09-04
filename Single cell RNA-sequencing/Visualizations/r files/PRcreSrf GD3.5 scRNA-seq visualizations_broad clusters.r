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

# Additional packages required for CellChat

library(CellChat);
library(patchwork);
options(stringsAsFactors=FALSE);
library(NMF);
library(ggalluvial);
library(ComplexHeatmap);
library(wordcloud);

# Record active package versions
packageVersion("Seurat");
packageVersion("ggplot2");        
packageVersion("dplyr");
packageVersion("viridis");
packageVersion("dittoSeq");
packageVersion("CellChat");
packageVersion("patchwork");
packageVersion("NMF");
packageVersion("ggalluvial");
packageVersion("ComplexHeatmap");
packageVersion("wordcloud");
packageVersion("scales");
packageVersion("scCustomize");

# Read in the processed RDS file from Sara Grimm, which holds a seurat object with filtered, integrated, annotated cells from 4 samples.
Adult_Srfflox_SrfKO <- readRDS("/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-FilteredIntegratedAnnotated.14mar2025.rds")
Adult_Srfflox_SrfKO
head(Adult_Srfflox_SrfKO)

DimPlot(Adult_Srfflox_SrfKO, group.by="Sample", pt.size=0.5, reduction='umap', shuffle=TRUE);

DimPlot(Adult_Srfflox_SrfKO, group.by="genotype", pt.size=0.5, reduction='umap', shuffle=TRUE);

DimPlot(Adult_Srfflox_SrfKO, group.by="Phase", pt.size=0.5, reduction='umap', shuffle=TRUE)

DimPlot(Adult_Srfflox_SrfKO, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE)

DimPlot(Adult_Srfflox_SrfKO, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

png(file = "PlotByCluster.png",
    width = 16, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend()

dev.off()

DimPlot(Adult_Srfflox_SrfKO, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype"
       )+NoLegend()

png(file = "PlotByCluster_SplitByGenotype.png",
    width = 29, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype"
        ) +NoLegend()

dev.off()

DimPlot(Adult_Srfflox_SrfKO, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )+NoLegend()

png(file = "PlotByBroadCellType_SplitByGenotype.png",
    width = 29, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
        )

dev.off()

png(file = "PlotByBroadCellType_SplitByGenotype_nolegend.png",
    width = 29, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
        )+NoLegend()

dev.off()

Idents(Adult_Srfflox_SrfKO) <- Adult_Srfflox_SrfKO$broad_cluster

FeaturePlot(Adult_Srfflox_SrfKO, features = c('Srf'), 
            pt.size=0.5,
            label=TRUE, 
            split.by="genotype",
            max.cutoff = 1.5,
            ) + theme(legend.position = "right")

png(file = "PlotSrf_SplitByGenotype_nomax.png",
    width = 29, height = 16, units = "cm", res = 600)

FeaturePlot(Adult_Srfflox_SrfKO, features = c('Srf'), 
            pt.size=0.5,
            label=TRUE, 
            split.by="genotype",
            ) + theme(legend.position = "right")

dev.off()

png(file = "PlotSrf_SplitByGenotype_max1-5.png",
    width = 29, height = 16, units = "cm", res = 600)

FeaturePlot(Adult_Srfflox_SrfKO, features = c('Srf'), 
            pt.size=0.5,
            label=TRUE, 
            split.by="genotype",
            max.cutoff = 1.5,
            ) + theme(legend.position = "right")

dev.off()

png(file = "PlotSrf_SplitByGenotype_max1-5_nolegend.png",
    width = 33, height = 16, units = "cm", res = 600)

FeaturePlot(Adult_Srfflox_SrfKO, features = c('Srf'), 
            pt.size=0.5,
            label=TRUE, 
            split.by="genotype",
            max.cutoff = 1.5,
            )

dev.off()

png(file = "PlotSrf_SplitByGenotype_max1-5_nolegend_nolabel.png",
    width = 33, height = 16, units = "cm", res = 600)

FeaturePlot(Adult_Srfflox_SrfKO, features = c('Srf'), 
            pt.size=0.5,
            label=FALSE, 
            split.by="genotype",
            max.cutoff = 1.5,
            )

dev.off()

FeaturePlot(Adult_Srfflox_SrfKO, features = c('Srf'), 
            pt.size=0.5,
            label=FALSE, 
            split.by="genotype",
            max.cutoff = 1.5,
            ) + theme(legend.position = "bottom", legend.text = element_text(angle = 0, size = 8))

png(file = "PlotSrf_SplitByGenotype_max1-5_nolegend_nolabel_turnscale.png",
    width = 33, height = 16, units = "cm", res = 600)

FeaturePlot(Adult_Srfflox_SrfKO, features = c('Srf'), 
            pt.size=0.5,
            label=FALSE, 
            split.by="genotype",
            max.cutoff = 1.5,
            ) + theme(legend.position = "bottom", legend.text = element_text(angle = 0, size = 8))

dev.off()

# Violin Plot Srf by broad cell type

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_tmp, features = c("Srf"), 
        group.by = "CellTypeBroad",
        split.by = "genotype",
        split.plot = TRUE,
        pt.size = 0.001,
        ncol = 1)
options(saved) # restore old settings

# Violin Plot Srf by broad cell type (smaller)

saved <- options(repr.plot.width = 4, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_tmp, features = c("Srf"), 
        group.by = "CellTypeBroad",
        split.by = "genotype",
        split.plot = TRUE,
        pt.size = 0.01,
        ncol = 1)
options(saved) # restore old settings

# Violin Plot Srf by broad cell type (turn)

saved <- options(repr.plot.width = 7, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_tmp, features = c("Srf"), 
        group.by = "CellTypeBroad",
        split.by = "genotype",
        split.plot = TRUE,
        pt.size = 0,
        ncol = 1) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
                          axis.text.y = element_text(angle = 45, vjust = 1, hjust=1),
                          axis.title.x = element_text(size = 0),
                          legend.text = element_text(angle = 45, size = 6),
                          legend.key.size = unit(1, "cm") ,
                          plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                          plot.title.position = "plot",
                          #aspect.ratio = 1/3
                          )

options(saved) # restore old settings

# Violin Plot Srf by broad cell type (turn), update colors

saved <- options(repr.plot.width = 7, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_tmp, features = c("Srf"), 
        group.by = "CellTypeBroad",
        split.by = "genotype",
        split.plot = TRUE,
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 1) + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1),
                          axis.text.y = element_text(angle = 45, vjust = 1, hjust=1),
                          axis.title.x = element_text(size = 0),
                          legend.text = element_text(angle = 45, size = 6),
                          legend.key.size = unit(1, "cm") ,
                          plot.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                          plot.title.position = "plot",
                          #aspect.ratio = 1/3
                          )

options(saved) # restore old settings

# Ridge Plot Srf by broad cell type

saved <- options(repr.plot.width = 12, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
RidgePlot(Adult_Srfflox_SrfKO_tmp, features = c("Srf"),
           #group.by = "genotype",
)
options(saved) # restore old settings

?VlnPlot

percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_tmp, features = "Srf", threshold = 0)
percent_stats

DimPlot(Adult_Srfflox_SrfKO, group.by="CellType", pt.size=0.5, reduction='umap', shuffle=TRUE)

DimPlot(Adult_Srfflox_SrfKO, group.by="CellSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE) 

DimPlot(Adult_Srfflox_SrfKO, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE) 

# Use DittoSeq to plot fraction of each cluster by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_tmp <- Adult_Srfflox_SrfKO
Idents(Adult_Srfflox_SrfKO_tmp) <- Adult_Srfflox_SrfKO_tmp$broad_cluster
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_tmp)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
# Use dittoBarPlot to generate the plot with some customization

saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "broad_cluster", 
             group.by = "genotype", 
             color.panel = my_color_palette, 
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings

# Change order
# Use DittoSeq to plot fraction of each cluster by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_tmp <- Adult_Srfflox_SrfKO
Idents(Adult_Srfflox_SrfKO_tmp) <- Adult_Srfflox_SrfKO_tmp$broad_cluster
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_tmp)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
# Use dittoBarPlot to generate the plot with some customization

saved <- options(repr.plot.width = 5, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "broad_cluster", 
             group.by = "genotype", 
             color.panel = my_color_palette, 
             var.labels.reorder = c(1, 2, 17, 18, 19, 13, 12, 15, 4, 7, 11, 14, 20, 16, 6, 10, 5, 3, 8, 9),
             main = "Cell Type Composition") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 3/1)

options(saved) # restore old settings

# Count cells by type to compare between genotypes
cell_table <- table(Adult_Srfflox_SrfKO_tmp$broad_cluster, Adult_Srfflox_SrfKO_tmp$genotype)
cell_table

# Use DittoSeq to plot fraction of each broad cell type by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_tmp <- Adult_Srfflox_SrfKO
Idents(Adult_Srfflox_SrfKO_tmp) <- Adult_Srfflox_SrfKO_tmp$CellTypeBroad
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_tmp)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
# Use dittoBarPlot to generate the plot with some customization

saved <- options(repr.plot.width = 8, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "CellTypeBroad", 
             group.by = "genotype",
             color.panel = my_color_palette, 
             main = "Broad Cell Type Composition") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 3/1)

options(saved) # restore old settings

# Use DittoSeq to plot fraction of each broad cell type by genotype (Turn)

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_tmp <- Adult_Srfflox_SrfKO
Idents(Adult_Srfflox_SrfKO_tmp) <- Adult_Srfflox_SrfKO_tmp$CellTypeBroad
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_tmp)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
# Use dittoBarPlot to generate the plot with some customization

saved <- options(repr.plot.width = 8, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "CellTypeBroad", 
             group.by = "genotype",
             color.panel = my_color_palette, 
             main = "") + theme(axis.text = element_text(size = 20, angle = 90),
                    axis.text.y = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    title = element_text(size = 20, hjust = 0.5), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    aspect.ratio = 3/1)

options(saved) # restore old settings


# Count cells by type to compare between genotypes
cell_table <- table(Adult_Srfflox_SrfKO_tmp$CellTypeBroad, Adult_Srfflox_SrfKO_tmp$genotype)
cell_table

# Use DittoSeq to plot fraction of each cell type by genotype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "genotype", 
             group.by = "CellTypeBroad", 
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 1/3)

options(saved) # restore old settings

# Use DittoSeq to plot fraction of each cell type by genotype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "genotype", 
             group.by = "CellType", 
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 1/3)

options(saved) # restore old settings

# Use DittoSeq to plot fraction of each cell type by genotype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "genotype", 
             group.by = "broad_cluster", 
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 1/3)

options(saved) # restore old settings

# Reorder, turn off label rotation, remove title
# Use DittoSeq to plot fraction of each cell type by genotype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_tmp, 
             var = "genotype", 
             group.by = "broad_cluster",
             x.reorder = c(1, 2, 17, 18, 19, 13, 12, 15, 4, 7, 11, 14, 20, 16, 6, 10, 5, 3, 8, 9),
             x.labels.rotate = FALSE,
             main = "") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 1/3)

options(saved) # restore old settings



head(Adult_Srfflox_SrfKO_tmp)
tail(Adult_Srfflox_SrfKO_tmp)

#Subset SrfKO only and Ctrl only
SrfKO_tmp <- subset(Adult_Srfflox_SrfKO_tmp, genotype == "SrfKO")

Ctrl_tmp <- subset(Adult_Srfflox_SrfKO_tmp, genotype == "Ctrl")

SrfKO_tmp
Ctrl_tmp


Idents(SrfKO_tmp) <- SrfKO_tmp$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = SrfKO_tmp, features = "Srf")
percent_stats

Idents(Ctrl_tmp) <- Ctrl_tmp$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Ctrl_tmp, features = "Srf")
percent_stats
