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
Adult_Srfflox_SrfKO_fibro <- readRDS("/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-Fibroblast.03feb2025.rds")
Adult_Srfflox_SrfKO_fibro
head(Adult_Srfflox_SrfKO_fibro)

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="broad_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

options(saved) # restore old settings

png(file = "Fibro_PlotBySubcluster.png",
    width = 16, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="seurat_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

dev.off()

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="seurat_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype"
       )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="seurat_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE, cols=c("#F8766D", "#B79F00", "#00BA38", "#00BFC4", "#619CFF", "#F564E3"), 
       )+NoLegend();

options(saved) # restore old settings

#adjust colors
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE, cols=c("#F8766D", #red
                           "#00BFC4", #aqua
                           "#619CFF", #blue
                           "#F564E3", #magenta
                           "#B79F00", #yellow
                           "#00BA38" #green
                          ))+NoLegend()

options(saved) # restore old settings

png(file = "Fibro_PlotBySubcluster_coloradjust.png",
    width = 16, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="seurat_clusters", pt.size=0.5, label.size = 5, reduction='umap', shuffle=TRUE,
        label=TRUE, cols=c("#F8766D", #red
                           "#00BFC4", #aqua
                           "#619CFF", #blue
                           "#F564E3", #magenta
                           "#B79F00", #yellow
                           "#00BA38" #green
                          ))+NoLegend()

dev.off()

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="FibroblastSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="FibroblastSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE
       );

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="FibroblastSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE
       )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="genotype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE
       )+NoLegend();

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="FibroblastSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="FibroblastSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )+NoLegend();

options(saved) # restore old settings

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$seurat_clusters

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="FibroblastSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype"
       )+NoLegend();

options(saved) # restore old settings

table(Adult_Srfflox_SrfKO_fibro$FibroblastSubtype,Adult_Srfflox_SrfKO_fibro$Sample);

table(Adult_Srfflox_SrfKO_fibro$FibroblastSubtype,Adult_Srfflox_SrfKO_fibro$genotype);

saved <- options(repr.plot.width = 8.5, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="Phase", pt.size=0.5, reduction='umap', shuffle=TRUE)

options(saved) # restore old settings

saved <- options(repr.plot.width = 15, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_fibro, group.by="Phase", split.by="genotype", pt.size=0.5, reduction='umap', shuffle=TRUE)

options(saved) # restore old settings

table(Adult_Srfflox_SrfKO_fibro$Phase,Adult_Srfflox_SrfKO_fibro$genotype);

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Pdgfra'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Mfap5'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hand2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Nr2f2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hoxa11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hoxa10'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")


FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Pgr'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Esr1'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Clec3b'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Col14a1'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Vim'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Mmp3'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

# Define marker genes

genes <- c("Pdgfra","Mfap5","Vim","Hoxa10","Pgr","Esr1", #pan fibroblast
           "Hand2","Hoxa11", #S1
           "Spon2","Aspg","Dpep1","Ngfr","Angptl7","Ifit1","Ifit3","H2-Q7", #F1
           "Cxcl14","Cdh11","Wt1","Rgs2","Smoc2","Wnt4","Aldh1a2","Bmp7", #F2
           "Clec3b","Col14a1","Fap","Cd55","Cxcl16","Mmp3","Efemp1","Vit", #F3
           "Top2a","Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cdc20","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67","Tpx2", # high in all of G2M
           "Hells","Ung","Cdc6","Dtl" # S phase
           );

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$seurat_clusters

saved <- options(repr.plot.width = 12, repr.plot.height = 3.5) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_fibro, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

# Define marker genes

genes <- c("Pdgfra","Mfap5","Vim","Hoxa10","Pgr","Esr1", #pan fibroblast
           "Hand2","Hoxa11", #S1
           "Spon2","Aspg","Dpep1","Ngfr","Angptl7","Ifit1","Ifit3","H2-Q7", #F1
           "Cxcl14","Cdh11","Wt1","Rgs2","Smoc2","Wnt4","Aldh1a2","Bmp7", #F2
           "Clec3b","Col14a1","Fap","Cd55","Cxcl16","Mmp3","Efemp1","Vit", #F3
           "Top2a","Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cdc20","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67","Tpx2", # high in all of G2M
           "Hells","Ung","Cdc6","Dtl" # S phase
           );

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_fibro, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

# Define marker genes

genes <- c("Pdgfra","Mfap5","Vim","Hoxa10","Pgr","Esr1", #pan fibroblast
           "Hand2","Hoxa11", #S1
           "Spon2","Aspg","Dpep1","Ngfr","Angptl7","Ifit1","Ifit3","H2-Q7", #F1
           "Cxcl14","Cdh11","Wt1","Rgs2","Smoc2","Wnt4","Aldh1a2","Bmp7", #F2
           "Clec3b","Col14a1","Fap","Cd55","Cxcl16","Mmp3","Efemp1","Vit", #F3
           "Top2a","Aurkb","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cdc20","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67","Tpx2", # high in all of G2M
           "Hells","Ung","Cdc6","Dtl" # S phase
           );

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 15, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_fibro, split.by="genotype", cols=c("#6b00ff","#6b00ff"), features=genes) + RotatedAxis() + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

# Define marker genes (reduce list)

genes <- c("Pdgfra", "Vim","Hoxa10","Pgr","Esr1", #pan fibroblast
           "Hand2","Hoxa11", #S1
           "Spon2","Aspg","Ngfr","Angptl7","H2-Q7", #F1
           "Cdh11","Wt1","Rgs2","Aldh1a2","Bmp7", #F2
           "Clec3b","Fap","Cd55","Cxcl16","Mmp3", #F3
           "Top2a","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67", # high in all of G2M
           "Hells","Dtl" # S phase
           );

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$seurat_clusters

saved <- options(repr.plot.width = 9, repr.plot.height = 3.5) # Make the plots bigger from here on out, save old options

DotPlot(Adult_Srfflox_SrfKO_fibro, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Pdgfra'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Vim'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hoxa10'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Pgr'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Esr1'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hand2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")


FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hoxa11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Spon2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Aspg'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Ngfr'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Angptl7'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('H2-Q7'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cdh11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Wt1'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Rgs2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Aldh1a2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Bmp7'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Clec3b'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Fap'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cd55'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cxcl16'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Mmp3'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Top2a'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Kif11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Ccnb2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cenpf'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Birc5'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Mki67'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hells'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Dtl'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings
      

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Pdgfra'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Vim'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hoxa10'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Pgr'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Esr1'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hand2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hoxa11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Spon2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Aspg'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Ngfr'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Angptl7'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('H2-Q7'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cdh11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Wt1'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Rgs2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Aldh1a2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Bmp7'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Clec3b'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Fap'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cd55'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cxcl16'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Mmp3'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Top2a'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Kif11'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Ccnb2'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Cenpf'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Birc5'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Mki67'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Hells'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = c('Dtl'), 
            pt.size=0.5,
            label=FALSE,
            repel = TRUE,
            )

options(saved) # restore old settings

genes <- c("Pdgfra", "Vim","Hoxa10","Pgr","Esr1", #pan fibroblast
           "Hand2","Hoxa11", #S1
           "Spon2","Aspg","Ngfr","Angptl7","H2-Q7", #F1
           "Cdh11","Wt1","Rgs2","Aldh1a2","Bmp7", #F2
           "Clec3b","Fap","Cd55","Cxcl16","Mmp3", #F3
           "Top2a","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67", # high in all of G2M
           "Hells","Dtl" # S phase
           );

saved <- options(repr.plot.width = 36, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = genes, 
            pt.size=0.5,
            label=FALSE,
            ncol=6
            )

options(saved) # restore old settings

genes <- c("Pdgfra", "Vim","Hoxa10","Pgr","Esr1", #pan fibroblast
           "Hand2","Hoxa11", #S1
           "Spon2","Ngfr","Angptl7", #F1
           "Rgs2","Aldh1a2","Bmp7", #F2
           "Clec3b","Cd55","Mmp3", #F3
           "Top2a","Kif11",  # G2M genes that are G2>M
           "Ccnb2","Cenpf",  # G2M genes that are M>G2
           "Birc5","Mki67", # high in all of G2M
           "Hells","Dtl" # S phase
           );

saved <- options(repr.plot.width = 36, repr.plot.height = 24) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = genes, 
            pt.size=0.5,
            label=FALSE,
            ncol=6
            )

options(saved) # restore old settings

genes <- c("Pdgfra", "Vim","Hoxa10","Pgr","Esr1", #pan fibroblast
           "Hand2","Hoxa11", #S1
           "Spon2","Ngfr","Angptl7", #F1
           "Rgs2","Bmp7", #F2
           "Clec3b","Mmp3", #F3
           "Top2a",  # G2M genes that are G2>M
           "Ccnb2",  # G2M genes that are M>G2
           "Mki67", # high in all of G2M
           "Hells" # S phase
           );

saved <- options(repr.plot.width = 36, repr.plot.height = 18) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = genes, 
            pt.size=0.5,
            label=FALSE,
            ncol=6
            )

options(saved) # restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 90) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_fibro, features = genes, 
            pt.size=0.5,
            label=FALSE,
            split.by="genotype",
            )

options(saved) # restore old settings


# Use DittoSeq to plot fraction of each cluster by genotype


saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_fibro, 
             var = "FibroblastSubtype", 
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

dittoBarPlot(Adult_Srfflox_SrfKO_fibro, 
             var = "FibroblastSubtype", 
             group.by = "genotype",
             color.panel = c("#F8766D", #red
                            "#B79F00", #yellow
                            "#00BA38", #green
                            "#F564E3", #magenta
                            "#619CFF"), #blue        #"#00BFC4", #aqua
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings   

# Use DittoSeq to plot fraction of each cluster by genotype (matching UMAP colors) - turn


saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_fibro, 
             var = "FibroblastSubtype", 
             group.by = "genotype",
             color.panel = c("#F8766D", #red
                            "#B79F00", #yellow
                            "#00BA38", #green
                            "#F564E3", #magenta
                            "#619CFF"), #blue        #"#00BFC4", #aqua
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20, angle = 90),
                    axis.text.y = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    title = element_text(size = 20, hjust = 0.5), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    aspect.ratio = 3/1)

#Subset proliferating fibroblasts

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

Adult_Srfflox_SrfKO_fibro_prolif <- subset(Adult_Srfflox_SrfKO_fibro, idents=c("ProliferatingFibroblastM","ProliferatingFibroblastSG2"));
Adult_Srfflox_SrfKO_fibro_prolif

table(Adult_Srfflox_SrfKO_fibro_prolif$FibroblastSubtype,Adult_Srfflox_SrfKO_fibro_prolif$genotype);

# Use DittoSeq to plot fraction of each cluster by genotype

saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_fibro_prolif, 
             var = "FibroblastSubtype", 
             group.by = "genotype",
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings      

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$seurat_clusters
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_fibro)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
my_color_palette

# Use DittoSeq to plot fraction of each cluster by genotype - alternate colors


saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_fibro_prolif, 
             var = "FibroblastSubtype", 
             group.by = "genotype",
             color.panel = c("#00BFC4", "#00BA38"),
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings      


# Use DittoSeq to plot fraction of each cluster by genotype - alternate colors


saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_fibro_prolif, 
             var = "FibroblastSubtype", 
             group.by = "genotype",
             color.panel = c("#F564E3", "#619CFF"),
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings      


# Use DittoSeq to plot fraction of each cluster by genotype - alternate colors (Turn)


saved <- options(repr.plot.width = 8, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_fibro_prolif, 
             var = "FibroblastSubtype", 
             group.by = "genotype",
             color.panel = c("#F564E3", "#619CFF"),
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20, angle = 90),
                    axis.text.y = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    axis.text.x = element_text(size = 20, angle = 90, vjust = 0.5, hjust=0.5),
                    title = element_text(size = 20, hjust = 0.5), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    aspect.ratio = 3/1)

options(saved) # restore old settings      

# Violin Plot inner stromal DEGs

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 24) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Zbp1", "Ifit3", "Ifi47", "Irf7", "Lyz2", "Oas1a", "Cxcl10", "Ifit1", "Mif", "Stat1", "Stat2", #inflammation
"Tgfbi", "Tgfb3", #Tgfb
"Mmp11", "Mmp2", #Mmps
"Col6a5", "Col6a4", "Col8a2", "Col7a1", "Col6a6", "Col25a1", "Col26a1", "Col1a1", "Col6a1", "Col6a2", "Col1a2", "Col6a3", #Colagens
"Sprr2f", "A2m", "Wnt4", "Wnt5a", "Ramp3", "Cxcl14", "Dio2", #Estrogen and senescence
"Mcm5", "Cdkn1c", "Hells", "Mcm3", "Mcm2", "Mcm4", "Mcm6", "Pcna" #proliferation
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot inner stromal DEGs (upregulated)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 48) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Zbp1", "Ifit3", "Ifi47", "Irf7", "Lyz2", "Oas1a", "Cxcl10", "Ifit1", "Mif", "Stat1", "Stat2", #inflammation
"Tgfbi", "Tgfb3", #Tgfb
"Mmp11", "Mmp2", #Mmps
"Col6a5", "Col6a4", "Col8a2", "Col7a1", "Col6a6", "Col25a1", "Col26a1", "Col1a1", "Col6a1", "Col6a2", "Col1a2", "Col6a3", #Colagens
"Sprr2f", "A2m", "Wnt4", "Wnt5a", "Ramp3", "Cxcl14", "Dio2", #Estrogen and senescence
"Mcm5", "Cdkn1c", "Hells", "Mcm3", "Mcm2", "Mcm4", "Mcm6", "Pcna" #proliferation
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot inner stromal DEGs (upregulated, reduce list)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 20) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1", "Stat2", #inflammation
"Tgfbi", #Tgfb
"Mmp11", "Mmp2", #Mmps
"Col6a4", "Col7a1", "Col6a6", "Col26a1", "Col1a1", "Col6a1", "Col6a2", "Col1a2", #Colagens
"A2m", "Wnt4", "Wnt5a", "Ramp3", "Cxcl14", "Dio2", #Estrogen and senescence
"Mcm5", "Cdkn1c", "Hells", "Mcm3", "Mcm2", "Mcm4", "Mcm6", "Pcna" #proliferation
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot inner stromal DEGs (downregulated)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 20) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1", "Stat2", #inflammation
"Tgfbi", #Tgfb
"Mmp11", "Mmp2", #Mmps
"Col6a4", "Col7a1", "Col6a6", "Col26a1", "Col1a1", "Col6a1", "Col6a2", "Col1a2", #Colagens
"A2m", "Wnt4", "Wnt5a", "Ramp3", "Cxcl14", "Dio2", #Estrogen and senescence
"Mcm5", "Cdkn1c", "Hells", "Mcm3", "Mcm2", "Mcm4", "Mcm6", "Pcna" #proliferation
"Egr1", "Fgf13", "Hgf", "Egfr", "Fos", "Myc", "Junb", "Stat3", #growth response
"Fst", "Scara5", "Zbtb16", "Fkbp5", "Klf15", #Progesterone
"Cenpf", "Mki67", #proliferation
"Tagln", "Actg1", "Myl9", "Acta2", "Vim", "Vcl", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot inner stromal DEGs (downregulated, reduce list)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 10) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Egr1", "Hgf", "Egfr", "Fos", "Myc", "Junb", "Stat3", #growth response/decidualization
"Fst", "Scara5", "Zbtb16", #Progesterone/decidualization
"Actg1", "Vim", "Vcl", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot inner stromal DEG selections

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 17.5) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1", "Stat2", #inflammation
"Col1a1", "Col1a2", "Tgfbi", "Mmp2", "Dio2", "A2m", "Cxcl14", #Cell stress and fibrosis
"Wnt4", "Wnt5a", "Ramp3", #Estrogen
"Egr1", "Hgf", "Egfr", "Fos", "Myc", "Junb", "Stat3", #growth response/decidualization
"Fst", "Scara5", "Zbtb16", #Progesterone/decidualization
"Actg1", "Vim", "Vcl", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot inner stromal DEG selections (reduced)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 15) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1", #inflammation
"Col1a1", "Col1a2", "Tgfbi", "Mmp2", "Dio2", "A2m", "Cxcl14", #Cell stress and fibrosis
"Egr1", "Hgf", "Egfr", "Fos", "Junb", "Stat3","Fst", "Scara5", "Zbtb16", #growth response/decidualization
"Actg1", "Vim", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot inner stromal DEG selections (reduced), remove some labels

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 12) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1", #inflammation
"Col1a1", "Col1a2", "Tgfbi", "Mmp2", "Dio2", "A2m", "Cxcl14", #Cell stress and fibrosis
"Egr1", "Hgf", "Egfr", "Fos", "Junb", "Stat3","Fst", "Scara5", "Zbtb16", #growth response/decidualization
"Actg1", "Vim", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings


# Violin Plot inner stromal DEG selections (reduced), remove some labels, change layout

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1", #inflammation
"Col1a1", "Col1a2", "Tgfbi", "Mmp2", "Dio2", "A2m", "Cxcl14", #Cell stress and fibrosis
"Egr1", "Hgf", "Egfr", "Fos", "Junb", "Stat3","Fst", "Scara5", "Zbtb16", #growth response/decidualization
"Actg1", "Vim", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings


# Violin Plot inner stromal DEG selections (reduced), remove some labels, change layout and drop a few

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Mmp2", "Dio2", "A2m", "Cxcl14", #Cell stress and fibrosis
"Isg15", "Mif", "Stat1", #inflammation
"Actg1", "Vim", "Actb", #cytoskeleton
"Egr1", "Hgf", "Egfr", "Fst", "Scara5", "Zbtb16" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings


# Violin Plot inner stromal DEG selections (reduced), remove some labels, change layout and drop a few

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Col6a1", "Col6a2", "Mmp2", "Dio2", #Cell stress and fibrosis
"Mif", "Isg15", "Stat1", #inflammation
"Actg1", "Vim", "Actb", #cytoskeleton
"Egr1", "Hgf", "Egfr", "Fst", "Zbtb16", "Scara5" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings


# Violin Plot inner stromal DEG selections (reduced), remove some labels, change layout and drop a few

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Col6a1", "Col6a2", "Mmp2", "Dio2", #Cell stress and fibrosis
"Mif", "Isg15", "Stat1", #inflammation
"Vim", "Actb", "Actg1",  #cytoskeleton
"Egr1", "Hgf", "Egfr", "Fst", "Zbtb16", "Scara5" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 6) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings



# Violin Plot inner stromal DEG selections (updated)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Col6a1", "Tgfbi", "Tgfb3", "Mmp2", "Dio2", #Cell stress and fibrosis
"Mif", "Isg15", "Stat1", #inflammation
"Vim", "Actb", "Actg1", "Vcl",  #cytoskeleton
"Egr1", "Fos", "Egfr", "Sgk1", "Fst", "Zbtb16", "Scara5" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings



# Violin Plot inner stromal DEG selections (updated)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Col6a1", "Tgfbi", "Tgfb3", "Mmp2", "Dio2", #Cell stress and fibrosis
"Mif", "Isg15", "Stat1", #inflammation
"Vim", "Actb", "Actg1", "Vcl",  #cytoskeleton
"Egr1", "Fos", "Egfr", "Hgf", "Fst", "Zbtb16", "Scara5" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings



# Violin Plot inner stromal DEG selections (updated)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Col6a1", "Tgfbi", "Tgfb3", "Mmp2", "Thbs1", #Cell stress and fibrosis
"Mif", "Isg15", "Stat1", #inflammation
"Vim", "Actb", "Actg1", "Vcl",  #cytoskeleton
"Egr1", "Fos", "Egfr", "Fst", "Zbtb16", "Scara5", "Dio2" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings



# Violin Plot inner stromal DEG selections (updated)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Col6a1", "Tgfbi", "Tgfb3", "Mmp2", "Snai1", #Cell stress and fibrosis
"Mif", "Isg15", "Stat1", #inflammation
"Vim", "Actb", "Actg1", "Vcl",  #cytoskeleton
"Egr1", "Fos", "Egfr", "Fst", "Zbtb16", "Scara5", "Dio2" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings



# Violin Plot inner stromal DEG selections (updated)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 14, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Col6a1", "Tgfbi", "Tgfb3", "Mmp2", "Snai1", #Cell stress and fibrosis
"Mif", "Isg15", "Stat1", #inflammation
"Vim", "Actb", "Actg1", "Vcl",  #cytoskeleton
"Egr1", "Fos", "Egfr", "Fst", "Zbtb16", "Scara5", "Clu" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 7) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings



# Violin Plot inner stromal DEG selections (reduced) - plot by row

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 2) # Make the plots bigger from here on out, save old options

VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1" #inflammation
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 9)

VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Col1a1", "Col1a2", "Tgfbi", "Mmp2", "Dio2", "A2m", "Cxcl14" #Cell stress and fibrosis
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 9)

VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Egr1", "Hgf", "Egfr", "Fos", "Junb", "Stat3","Fst", "Scara5", "Zbtb16" #growth response/decidualization
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 9)

VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Actg1", "Vim", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 9)

options(saved) # restore old settings

  # Violin Plot inner stromal DEG selections (reduced) - stack

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 5, repr.plot.height = 10) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Isg15", "Mif", "Stat1", #inflammation
"Col1a1", "Col1a2", "Tgfbi", "Mmp2", "Dio2", "A2m", "Cxcl14", #Cell stress and fibrosis
"Egr1", "Hgf", "Egfr", "Fos", "Junb", "Stat3","Fst", "Scara5", "Zbtb16", #growth response/decidualization
"Actg1", "Vim", "Actb" #cytoskeleton
), 
        idents = "FibroblastF2", 
        group.by = "genotype",
        pt.size = 0,
        stack = TRUE,
        flip = TRUE,

)
options(saved) # restore old settings



#Subset inner stromal fibroblasts

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

Adult_Srfflox_SrfKO_fibro_strom <- subset(Adult_Srfflox_SrfKO_fibro, idents=c("FibroblastF2"));
Adult_Srfflox_SrfKO_fibro_strom

# Single cell heatmap of feature expression
Idents(Adult_Srfflox_SrfKO_fibro_strom) <- Adult_Srfflox_SrfKO_fibro_strom$genotype

saved <- options(repr.plot.width = 8, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

DoHeatmap(subset(Adult_Srfflox_SrfKO_fibro_strom, downsample = 100), features = c(
"Isg15", "Mif", "Stat1", "Stat2", #inflammation
"Tgfbi", #Tgfb
"Mmp11", "Mmp2", #Mmps
"Col6a4", "Col7a1", "Col6a6", "Col26a1", "Col1a1", "Col6a1", "Col6a2", "Col1a2", #Colagens
"A2m", "Wnt4", "Wnt5a", "Ramp3", "Cxcl14", "Dio2", #Estrogen and senescence
"Mcm5", "Cdkn1c", "Hells", "Mcm3", "Mcm2", "Mcm4", "Mcm6", "Pcna", #proliferation
"Egr1", "Hgf", "Egfr", "Fos", "Myc", "Junb", "Stat3", #growth response/decidualization
"Fst", "Scara5", "Zbtb16", #Progesterone/decidualization
"Actg1", "Vim", "Vcl", "Actb" #cytoskeleton
), size = 3)

options(saved) # restore old settings

# Violin Plot outer basal DEG selections

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 22.5) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Mif", "Ifi35", "Cxcl16", "Ifitm1", "Apoe", "Lyz2", "Irf1", "Cxcl12", "C2", "Isg15", "Cxcl10", #inflammation
"Col12a1", "Col11a1", "Col4a5", "Col5a2", "Mmp2", "Mmp11", #Cell stress and fibrosis
"Egr1", "Zbtb16",   #growth response/decidualization
"Vegfa", "Vcam1", "Vegfb", #vascular
"Greb1", "Ramp3", #E2
"Acta2", "Tagln", "Cnn2", "Actb", "Actg1", "Myh9", "Myl6", "Myl9", "Vcl", "Tpm2", "Tpm4", "Vit" #cytoskeleton
), 
        idents = "FibroblastF3", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot outer basal DEG selections (reduced)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 17.5) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Mif", "Ifi35", "Cxcl16", "Ifitm1", "Irf1", "Cxcl12", "C2", #inflammation
"Col12a1", "Col11a1", "Col4a5", "Col5a2", "Mmp2", "Mmp11", #Cell stress and fibrosis
"Egr1", "Zbtb16",   #growth response/decidualization
"Vegfa", "Vcam1", "Vegfb", #vascular
"Acta2", "Cnn2", "Actb", "Actg1", "Myh9", "Myl6", "Myl9", "Vcl", "Tpm2", "Tpm4" #cytoskeleton
), 
        idents = "FibroblastF3", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot outer basal DEG selections (reduced-2)

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 10) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Mif", "Ifitm1", "Irf1", #inflammation
"Mmp2", #Cell stress and fibrosis
"Egr1", "Zbtb16",   #growth response/decidualization
"Acta2", "Cnn2", "Actb", "Actg1", "Myh9", "Myl6", "Myl9", "Vcl", "Tpm2", "Tpm4" #cytoskeleton
), 
        idents = "FibroblastF3", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4)
options(saved) # restore old settings

# Violin Plot outer basal DEG selections (reduced-2), remove some labels

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Mif", "Ifitm1", "Irf1", #inflammation
"Mmp2", #Cell stress and fibrosis
"Egr1", "Zbtb16",   #growth response/decidualization
"Acta2", "Cnn2", "Actb", "Actg1", "Myh9", "Myl6", "Myl9", "Vcl", "Tpm2", "Tpm4" #cytoskeleton
), 
        idents = "FibroblastF3", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 4) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings

# Violin Plot outer basal DEG selections (reduced-2), remove some labels

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Mif", "Ifitm1", "Irf1", #inflammation
"Mmp2", #Cell stress and fibrosis
"Egr1", "Zbtb16",   #growth response/decidualization
"Acta2", "Cnn2", "Actb", "Actg1", "Myh9", "Myl6", "Myl9", "Vcl", "Tpm2", "Tpm4" #cytoskeleton
), 
        idents = "FibroblastF3", 
        group.by = "genotype",
        pt.size = 0,
        ncol = 6) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings

# Violin Plot outer basal DEG selections (reduced-2), remove some labels, change colors

Idents(Adult_Srfflox_SrfKO_fibro) <- Adult_Srfflox_SrfKO_fibro$FibroblastSubtype

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
VlnPlot(Adult_Srfflox_SrfKO_fibro, features = c(
"Mif", "Ifitm1", "Irf1", #inflammation
"Mmp2", #Cell stress and fibrosis
"Egr1", "Zbtb16",   #growth response/decidualization
"Acta2", "Cnn2", "Actb", "Actg1", "Myh9", "Myl6", "Myl9", "Vcl", "Tpm2", "Tpm4" #cytoskeleton
), 
        idents = "FibroblastF3", 
        group.by = "genotype",
        cols = c("dimgrey", "firebrick2"),
        pt.size = 0,
        ncol = 6) & theme(axis.title=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
options(saved) # restore old settings
