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

# Read in the processed RDS file from Sara Grimm, which holds a seurat object with filtered, integrated, annotated cells from 4 samples.
# We are using the broad clusters for these cell types
Adult_Srfflox_SrfKO <- readRDS("/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-FilteredIntegratedAnnotated.14mar2025.rds")
Adult_Srfflox_SrfKO
head(Adult_Srfflox_SrfKO)

DimPlot(Adult_Srfflox_SrfKO, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype"
       )+NoLegend()

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

Adult_Srfflox_SrfKO_tmp <- Adult_Srfflox_SrfKO
Idents(Adult_Srfflox_SrfKO_tmp) <- Adult_Srfflox_SrfKO_tmp$broad_cluster

DimPlot(Adult_Srfflox_SrfKO_tmp, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )

options(saved) # restore old settings

Idents(Adult_Srfflox_SrfKO_tmp) <- Adult_Srfflox_SrfKO_tmp$broad_cluster

SM_PV_cells <- WhichCells(Adult_Srfflox_SrfKO_tmp, idents = c(4,11,14))

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
DimPlot(object = Adult_Srfflox_SrfKO_tmp, cells = SM_PV_cells, reduction = "umap", label = T)
options(saved) # restore old settings

#Subset clusters 4, 11, 14 only

Idents(Adult_Srfflox_SrfKO_tmp) <- Adult_Srfflox_SrfKO_tmp$broad_cluster

Adult_Srfflox_SrfKO_SM_PV <- subset(Adult_Srfflox_SrfKO_tmp, idents=c(4,11,14));
Adult_Srfflox_SrfKO_SM_PV


table(Adult_Srfflox_SrfKO_SM_PV$broad_cluster,Adult_Srfflox_SrfKO_SM_PV$Sample);

table(Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad,Adult_Srfflox_SrfKO_SM_PV$genotype);

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_SM_PV, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_SM_PV, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )

options(saved) # restore old settings

saved <- options(repr.plot.width = 15, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_SM_PV, group.by="broad_cluster", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
            )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_SM_PV, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )

options(saved) # restore old settings

saved <- options(repr.plot.width = 15, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_SM_PV, group.by="CellTypeBroad", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
            )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 15, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_SM_PV, group.by="Phase", split.by="genotype", pt.size=0.5, reduction='umap', shuffle=TRUE)

options(saved) # restore old settings

table(Adult_Srfflox_SrfKO_SM_PV$Phase,Adult_Srfflox_SrfKO_SM_PV$genotype);

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Acta2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Tagln'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Myh11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Cnn1'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Actg2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Pdgfrb'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Rgs5'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Mcam'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_SM_PV, features = c('Nes'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Acta2",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Tagln",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Myh11",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

SMC_marker_pos <- WhichCells(object = Adult_Srfflox_SrfKO_SM_PV, expression = Acta2 > 0 & Tagln > 0 & Myh11 > 0)
SMC_marker_pos

Adult_Srfflox_SrfKO_SM_PV_SMC_marker_pos <- subset(Adult_Srfflox_SrfKO_SM_PV, cells = SMC_marker_pos);
Adult_Srfflox_SrfKO_SM_PV_SMC_marker_pos

table(Adult_Srfflox_SrfKO_SM_PV_SMC_marker_pos$CellTypeBroad,Adult_Srfflox_SrfKO_SM_PV_SMC_marker_pos$genotype);

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Cnn1",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Actg2",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Pdgfrb",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Mcam",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Nes",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
percent_stats <- Percent_Expressing(seurat_object = Adult_Srfflox_SrfKO_SM_PV, 
                                    features = "Rgs5",
                                    group_by = "CellTypeBroad",
                                    split_by = "genotype",
                                   )
percent_stats

genes <- c("Acta2","Tagln","Myh11","Cnn1","Actg2","Pdgfrb", "Rgs5", "Mcam","Nes","Susd2");
DotPlot(Adult_Srfflox_SrfKO_SM_PV, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

genes <- c("Acta2","Tagln","Myh11","Cnn1","Actg2","Pdgfrb", "Rgs5", "Mcam","Nes","Susd2");
DotPlot(Adult_Srfflox_SrfKO_SM_PV, split.by="genotype", cols=c("forestgreen","mediumblue"), features=genes) + RotatedAxis() + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

# Tried to change order of plotting as below but it didn't work (https://www.r-bloggers.com/2007/10/reorder-factor-levels/)

# check
# levels(factor("CellTypeBroad"))
# levels(factor("genotype"))
#Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$CellTypeBroad
#levels(Adult_Srfflox_SrfKO_SM_PV)
#Idents(Adult_Srfflox_SrfKO_SM_PV) <- Adult_Srfflox_SrfKO_SM_PV$genotype
#levels(Adult_Srfflox_SrfKO_SM_PV)

#Dropped two of the rows from graph, not sure why
## To reorder the levels:
## note, if x is not a factor use levels(factor(x))
# "CellTypeBroad" = factor("CellTypeBroad",levels("CellTypeBroad")[c("SmoothMuscle","Perivascular",1:2)])
# "genotype" = factor("genotype",levels("genotype")[c("Ctrl","SrfKO",1:2)])

#levels(factor("CellTypeBroad", c("SmoothMuscle","Perivascular",1:2)))
#levels(factor("genotype", c("Ctrl","SrfKO",1:2)))

saved <- options(repr.plot.width = 12, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

genes <- c("Acta2","Tagln","Myh11","Cnn1","Actg2","Pdgfrb", "Rgs5", "Mcam","Nes");
DotPlot(Adult_Srfflox_SrfKO_SM_PV, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

saved <- options(repr.plot.width = 10, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

genes <- c("Acta2","Tagln","Myh11","Cnn1","Actg2","Pdgfrb", "Rgs5", "Mcam","Nes");
DotPlot(Adult_Srfflox_SrfKO_SM_PV, split.by="genotype", cols=c("#6b00ff","#6b00ff"), features=genes) + RotatedAxis() + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

saved <- options(repr.plot.width = 9, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

genes <- c("Acta2","Tagln","Myh11","Cnn1","Actg2","Pdgfrb", "Rgs5","Nes");
DotPlot(Adult_Srfflox_SrfKO_SM_PV, split.by="genotype", cols=c("#6b00ff","#6b00ff"), features=genes) + RotatedAxis() + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings
