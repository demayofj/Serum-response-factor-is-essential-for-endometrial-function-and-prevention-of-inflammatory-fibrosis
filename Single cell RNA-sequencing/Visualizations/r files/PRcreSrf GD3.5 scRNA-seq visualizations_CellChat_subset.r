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
library(presto)

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
packageVersion("presto");

#Read in the data from Sara Grimm
load(file = "/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-CellChatDB-merged.20feb2025.RData")

# Contains "CCobjMerged", a merged CellChat object containing CCobjCtrl and CCobjSrfKO CellChat objects that were first run separately.

CCobjMerged



ls()

object.list

alphaOrderCellTypes

CCobjMerged@options

CCobjMerged@meta

levels(CCobjCtrl@idents)
levels(CCobjSrfKO@idents)

#Try subsetting cell types of interest
#see https://rdrr.io/github/sqjin/CellChat/man/subsetCellChat.html

#First Ctrl
CCobjCtrl_sub <-subsetCellChat(
  CCobjCtrl,
  cells.use = NULL,
  idents.use = c('Epi-Luminal','Epi-Glnd','Fibro-F2','Fibro-F3',
                 'Imm-Macrophage','Imm-Neutrophil'),
  group.by = NULL,
  invert = FALSE,
  thresh = 0.05
)

#Also the SrfKO object

CCobjSrfKO_sub <-subsetCellChat(
  CCobjSrfKO,
  cells.use = NULL,
  idents.use = c('Epi-Luminal','Epi-Glnd','Fibro-F2','Fibro-F3',
                 'Imm-Macrophage','Imm-Neutrophil'),
  group.by = NULL,
  invert = FALSE,
  thresh = 0.05
)

levels(CCobjCtrl_sub@idents);
levels(CCobjSrfKO_sub@idents);

#Re-merge

object.list_sub <- list(Ctrl_sub=CCobjCtrl_sub, SrfKO_sub=CCobjSrfKO_sub);
CCobjMerged_sub <- mergeCellChat(object.list_sub, add.names=names(object.list_sub));

object.list_sub

CCobjMerged_sub@meta

## Part I: Predict general principles of cell-cell communication.

saved <- options(repr.plot.width = 5, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

compareInteractions(CCobjMerged_sub, show.legend=F, group=c(1,2))

options(saved) # restore old settings

saved <- options(repr.plot.width = 5, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

compareInteractions(CCobjMerged_sub, show.legend=F, group=c(1,2), measure="weight")

options(saved) # restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(CCobjMerged_sub, weight.scale = T)
netVisual_diffInteraction(CCobjMerged_sub, weight.scale = T, measure = "weight")

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netVisual_heatmap(CCobjMerged_sub)

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netVisual_heatmap(CCobjMerged_sub, measure="weight")

options(saved) # restore old settings

#adjust colors


my_color_palette = c("Epi-Glnd"="gold", "Epi-Luminal"="mediumorchid", "Fibro-F2"="#F8766D", "Fibro-F3"="#B79F00", "Imm-Macrophage"="seagreen2", "Imm-Neutrophil"="dodgerblue2")


saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netVisual_heatmap(CCobjMerged_sub, color.use = my_color_palette, measure="weight")

options(saved) # restore old settings


?netVisual_heatmap

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

weight.max <- getMaxWeight(object.list_sub, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_sub)) {
  netVisual_circle(object.list_sub[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list_sub)[i]))
}

options(saved) # restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

weight.max <- getMaxWeight(object.list_sub, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list_sub)) {
  netVisual_circle(object.list_sub[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Weight of interactions - ", names(object.list_sub)[i]))
}

options(saved) # restore old settings

## compare major sources and targets in 2D space

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

num.link <- sapply(object.list_sub, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)});
weight.MinMax <- c(min(num.link), max(num.link));
gg <- list();
for (i in 1:length(object.list_sub)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list_sub[[i]], title=names(object.list_sub)[i], weight.MinMax=weight.MinMax);
}
patchwork::wrap_plots(plot=gg);

options(saved) # restore old settings



gg1 <- rankNet(CCobjMerged_sub, mode="comparison", stacked=T, do.stat=TRUE,color.use=c("dimgrey", "firebrick2"),);
gg2 <- rankNet(CCobjMerged_sub, mode="comparison", stacked=F, do.stat=TRUE,color.use=c("dimgrey", "firebrick2"),);

saved <- options(repr.plot.width = 10, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1+gg2

options(saved) #restore old settings


?rankNet

#Reduce list by selecting signaling types
gg1 <- rankNet(CCobjMerged_sub, mode="comparison", stacked=T, do.stat=TRUE, color.use=c("dimgrey", "firebrick2"),
               signaling = c("EGF", "COMPLEMENT", "MMP", "MIF", "LIFR", "FGF", "Prostaglandin", "FN1", "COLLAGEN", "TNF", 
                             "CSF", "CXCL", "TGFb", "WNT", "CX3C", "CSF3", "RBP4", "BMP"))

saved <- options(repr.plot.width = 5, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

gg1

options(saved) #restore old settings


#flip
gg1 <- rankNet(CCobjMerged_sub, mode="comparison", stacked=T, do.stat=TRUE, color.use=c("dimgrey", "firebrick2"), do.flip = F,
               signaling = c("EGF", "COMPLEMENT", "MMP", "MIF", "LIFR", "FGF", "Prostaglandin", "FN1", "COLLAGEN", "TNF", 
                             "CSF", "CXCL", "TGFb", "WNT", "CX3C", "CSF3"))

saved <- options(repr.plot.width = 5, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

gg1

options(saved) #restore old settings


levels(CCobjSrfKO_sub@idents);

#1 'Epi-Glnd'
#2 'Epi-Luminal'
#3 'Fibro-F2'
#4 'Fibro-F3'
#5 'Imm-Macrophage'
#6 'Imm-Neutrophil'

## Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs.
## NOTE FROM SARA -- Most plots in part III require specifying the source & target cell types ... 
# I'm only making example views here, you'll want to modify the sources.use and targets.use parameters according to your interests.

saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties

netVisual_bubble(CCobjMerged_sub, sources.use=2, targets.use=c(5:6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


#Confirmed that it matches the full set - try one more just to be safe
saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties

netVisual_bubble(CCobjMerged_sub, sources.use=5, targets.use=c(3:4), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


## Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram.
## NOTE FROM SARA -- Most plots in part IV require specifying pathways and sometimes also the source & target cell types; you'll want to modify the signaling/sources.use/targets.use parameters according to your interests.

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("CXCL"); 
weight.max <- getMaxWeight(object.list_sub, slot.name=c("netP"), attribute=pathways.show);
par(mfrow=c(1,1), xpd=TRUE);
netVisual_aggregate(object.list_sub[[1]], signaling=pathways.show, layout="circle", edge.weight.max=weight.max[1], edge.width.max=10, signaling.name=paste(pathways.show, names(object.list_sub)[1]));

options(saved) #restore old settings

## Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram.
## NOTE FROM SARA -- Most plots in part IV require specifying pathways and sometimes also the source & target cell types; you'll want to modify the signaling/sources.use/targets.use parameters according to your interests.

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("MIF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

#adjust colors, remove labels by making them tiny

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

my_color_palette = c("Epi-Glnd"="gold2", "Epi-Luminal"="mediumorchid", "Fibro-F2"="#F8766D", "Fibro-F3"="darkgoldenrod4", "Imm-Macrophage"="seagreen2", "Imm-Neutrophil"="dodgerblue2")

pathways.show <- c("MIF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, color.use = my_color_palette,   vertex.label.cex = 0.000000001,
                      signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("CCL") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("COLLAGEN") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

#show only signal coming from macrophages and neutrophils

pathways.show <- c("COLLAGEN") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, sources.use = c("Imm-Macrophage", "Imm-Neutrophil"), layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings



saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("FN1") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("COMPLEMENT") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

#adjust colors

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

my_color_palette = c("Epi-Glnd"="gold", "Epi-Luminal"="mediumorchid", "Fibro-F2"="#F8766D", "Fibro-F3"="#B79F00", "Imm-Macrophage"="seagreen2", "Imm-Neutrophil"="dodgerblue2")

pathways.show <- c("COMPLEMENT") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, color.use = my_color_palette, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

#adjust colors, remove labels by making them tiny

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

my_color_palette = c("Epi-Glnd"="gold2", "Epi-Luminal"="mediumorchid", "Fibro-F2"="#F8766D", "Fibro-F3"="darkgoldenrod4", "Imm-Macrophage"="seagreen2", "Imm-Neutrophil"="dodgerblue2")

pathways.show <- c("COMPLEMENT") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, color.use = my_color_palette,   vertex.label.cex = 0.000000001,
                      signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("CSF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("EGF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

#adjust colors, remove labels by making them tiny

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

my_color_palette = c("Epi-Glnd"="gold2", "Epi-Luminal"="mediumorchid", "Fibro-F2"="#F8766D", "Fibro-F3"="darkgoldenrod4", "Imm-Macrophage"="seagreen2", "Imm-Neutrophil"="dodgerblue2")

pathways.show <- c("EGF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, color.use = my_color_palette,   vertex.label.cex = 0.000000001,
                      signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("FGF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("IL1") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings


saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

#Show only signal coming from macrophages

pathways.show <- c("IL1") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, sources.use = "Imm-Macrophage", signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings


saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings


#adjust colors, remove labels by making them tiny

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

my_color_palette = c("Epi-Glnd"="gold2", "Epi-Luminal"="mediumorchid", "Fibro-F2"="#F8766D", "Fibro-F3"="darkgoldenrod4", "Imm-Macrophage"="seagreen2", "Imm-Neutrophil"="dodgerblue2")

pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], 
                      edge.width.max = 10, color.use = my_color_palette,   vertex.label.cex = 0.000000001,
                      signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings


saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options

pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list_sub[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list_sub)[i]))
}

options(saved) #restore old settings


pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow=c(1,1), xpd=TRUE);
netVisual_aggregate(object.list_sub[[1]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list_sub)[1]));


netVisual_aggregate(object.list_sub[[2]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list_sub)[2]));

pathways.show <- c("MIF") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow=c(1,1), xpd=TRUE);
netVisual_aggregate(object.list_sub[[1]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list_sub)[1]));


netVisual_aggregate(object.list_sub[[2]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list_sub)[2]));

pathways.show <- c("COMPLEMENT") 
weight.max <- getMaxWeight(object.list_sub, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow=c(1,1), xpd=TRUE);
netVisual_aggregate(object.list_sub[[1]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list_sub)[1]));


netVisual_aggregate(object.list_sub[[2]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list_sub)[2]));

# All the signaling pathways showing significant communications can be accessed by
CCobjSrfKO_sub@netP$pathways

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CXCL") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CXCL <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CXCL <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("COMPLEMENT") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COMPLEMENT <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COMPLEMENT[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COMPLEMENT <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COMPLEMENT[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("COLLAGEN") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[9,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[10,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("MMP") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.MMP <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MMP[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("IL1") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.IL1 <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IL1[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.IL1 <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IL1[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("MIF") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.MIF <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MIF[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.MIF <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MIF[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CSF") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CSF <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CSF[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CCL") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[10,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("FN1") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.FN1 <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.FN1[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.FN1 <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.FN1[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CX3C") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CX3C <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CX3C[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("TGFb") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[3,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[4,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[5,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[6,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[7,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[8,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[9,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("EGF") 
netAnalysis_contribution(CCobjSrfKO_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGF <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGF <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGF <- extractEnrichedLR(CCobjSrfKO_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[3,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("EGF") 
netAnalysis_contribution(CCobjCtrl_sub, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGF <- extractEnrichedLR(CCobjCtrl_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjCtrl_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGF <- extractEnrichedLR(CCobjCtrl_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjCtrl_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGF <- extractEnrichedLR(CCobjCtrl_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[3,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjCtrl_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.EGF <- extractEnrichedLR(CCobjCtrl_sub, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.EGF[4,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjCtrl_sub, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
