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

#Read in the data
load(file = "/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-CellChatDB-merged.20feb2025.RData")

# Contains "CCobjMerged", a merged CellChat object containing CCobjCtrl and CCobjSrfKO CellChat objects that were first run separately.

CCobjMerged



ls()

object.list

alphaOrderCellTypes

CCobjMerged@options

CCobjMerged@meta

## Part I: Predict general principles of cell-cell communication.

saved <- options(repr.plot.width = 5, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

compareInteractions(CCobjMerged, show.legend=F, group=c(1,2))

options(saved) # restore old settings

saved <- options(repr.plot.width = 5, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

compareInteractions(CCobjMerged, show.legend=F, group=c(1,2), measure="weight")

options(saved) # restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(CCobjMerged, weight.scale = T)
netVisual_diffInteraction(CCobjMerged, weight.scale = T, measure = "weight")

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netVisual_heatmap(CCobjMerged)

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netVisual_heatmap(CCobjMerged, measure="weight")

options(saved) # restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

options(saved) # restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Weight of interactions - ", names(object.list)[i]))
}

options(saved) # restore old settings

## compare major sources and targets in 2D space

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)});
weight.MinMax <- c(min(num.link), max(num.link));
gg <- list();
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title=names(object.list)[i], weight.MinMax=weight.MinMax);
}
patchwork::wrap_plots(plot=gg);

options(saved) # restore old settings



alphaOrderCellTypes

## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="Epi-Luminal")

options(saved) # restore old settings


## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="Epi-Glnd")

options(saved) # restore old settings


## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="Fibro-F2")

options(saved) # restore old settings

## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="Fibro-F3")

options(saved) # restore old settings

## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="Imm-Macrophage")

options(saved) # restore old settings

## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="Imm-Neutrophil")

options(saved) # restore old settings

## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="Perivascular")

options(saved) # restore old settings

## plot differential signaling scatterplot for user-defined cell type

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netAnalysis_signalingChanges_scatter(CCobjMerged, idents.use="SmoothMuscle")

options(saved) # restore old settings

## Part II: Identify the conserved and context-specific signaling pathways.

CCobjMerged <- computeNetSimilarityPairwise(CCobjMerged, type="functional");
CCobjMerged <- netEmbedding(CCobjMerged, type="functional", umap.method="uwot");
CCobjMerged <- netClustering(CCobjMerged, type="functional");

netVisual_embeddingPairwise(CCobjMerged, type="functional", label.size=3.5)


saved <- options(repr.plot.width = 10, repr.plot.height = 15) # Make the plots bigger from here on out, save old options

netVisual_embeddingPairwiseZoomIn(CCobjMerged, type="functional", nCol=2);

options(saved) # restore old settings


#Warning, this one takes some time (30-60 min?)
CCobjMerged <- computeNetSimilarityPairwise(CCobjMerged, type="structural");
CCobjMerged <- netEmbedding(CCobjMerged, type="structural", umap.method="uwot");
CCobjMerged <- netClustering(CCobjMerged, type="structural");
netVisual_embeddingPairwise(CCobjMerged, type="structural", label.size=3.5)

netVisual_embeddingPairwiseZoomIn(CCobjMerged, type="structural", nCol=3);

saved <- options(repr.plot.width = 8, repr.plot.height = 24) # Make the plots bigger from here on out, save old options

rankSimilarity(CCobjMerged, type="functional");  # only evaluates signaling pathways identified in both datasets

options(saved) # restore old settings


rankSimilarity(CCobjMerged, type="structural");  # only evaluates signaling pathways identified in both datasets

gg1 <- rankNet(CCobjMerged, mode="comparison", stacked=T, do.stat=TRUE);
gg2 <- rankNet(CCobjMerged, mode="comparison", stacked=F, do.stat=TRUE);

saved <- options(repr.plot.width = 10, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1+gg2

options(saved) #restore old settings


saved <- options(repr.plot.width = 30, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

i=1;
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways);
png(file = "cellchat_outputs/merged.signalingRole-outgoing.heatmap-view.png",
    width = 16, height = 16, units = "cm", res = 600)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern="outgoing", signaling=pathway.union, title=names(object.list)[i], width=8, height=20);
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern="outgoing", signaling=pathway.union, title=names(object.list)[i+1], width=8, height=20);
draw(ht1+ht2, ht_gap=unit(0.5, "cm"));
dev.off();

options(saved) #restore old settings


saved <- options(repr.plot.width = 30, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

i=1;
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways);
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern="outgoing", signaling=pathway.union, title=names(object.list)[i], width=8, height=20);
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern="outgoing", signaling=pathway.union, title=names(object.list)[i+1], width=8, height=20);
draw(ht1+ht2, ht_gap=unit(0.5, "cm"));

options(saved) #restore old settings


alphaOrderCellTypes

#1 'EndoLymph'
#2 'Endothelial'
#3 'Epi-Glnd'
#4 'Epi-Glnd_Prolif'
#5 'Epi-LowQual'
#6 'Epi-Luminal'
#7 'Fibro-F2'
#8 'Fibro-F3'
#9 'Fibro-ProlifG2M'
#10 'Fibro-ProlifSG2'
#11 'Fibro-Unk'
#12 'Imm-Bcell'
#13 'Imm-cDC2_mregDC'
#14 'Imm-DC1'
#15 'Imm-gdT'
#16 'Imm-ILC2'
#17 'Imm-Macrophage'
#18 'Imm-Monocyte'
#19 'Imm-Neutrophil'
#20 'Imm-NK'
#21 'Imm-NK_Prolif'
#22 'Imm-pDC'
#23 'Imm-Plasma'
#24 'Imm-T_NKT'
#25 'Mesothelial'
#26 'Perivascular'
#27 'SmoothMuscle'

## Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs.
## NOTE FROM SARA -- Most plots in part III require specifying the source & target cell types ... 
# I'm only making example views here, you'll want to modify the sources.use and targets.use parameters according to your interests.

saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
#7 'Fibro-F2'
#17 'Imm-Macrophage'
#18 'Imm-Monocyte'
#19 'Imm-Neutrophil'
#20 'Imm-NK'
#21 'Imm-NK_Prolif'
#22 'Imm-pDC'

netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(17:22), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17:19), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 5, repr.plot.height = 25) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:10), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 5, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:8), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=19, targets.use=c(7:10), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(17:19), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 5, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(17,19), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=4, targets.use=c(17:19), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=18, targets.use=c(7:10), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=18, targets.use=c(3,4,6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(3,4,6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=19, targets.use=c(3,4,6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 50) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(3,4,6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 50) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=10, targets.use=c(3,4,6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 50) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=9, targets.use=c(3,4,6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 50) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(7:10), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 50) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:10), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 5, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:8), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 10, repr.plot.height = 50) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=4, targets.use=c(7:10), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


saved <- options(repr.plot.width = 5, repr.plot.height = 40) # Make the plots bigger from here on out, save old options

## identify the up- and down-regulated signaling ligand-receptor pairs based on communication probabilties
alphaOrderCellTypes <- sort(unique(CCobjCtrl@idents));
netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(3,6), comparison=c(1,2), angle.x=45);

options(saved) #restore old settings


#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 

saved <- options(repr.plot.width = 10, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(17:19), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(17:19), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17:19), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17:19), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 10, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

#limit to complement

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("COMPLEMENT"))

saved <- options(repr.plot.width = 5, repr.plot.height = 4.9) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=F);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

#limit to complement

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("COMPLEMENT"))

saved <- options(repr.plot.width = 5, repr.plot.height = 4.9) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=c(3,6), targets.use=c(17,19), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=F);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

#limit to complement

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("COMPLEMENT"))

saved <- options(repr.plot.width = 6, repr.plot.height = 4.9) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=c(4,6,8,17,19), targets.use=c(17,19), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=F);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

#limit to complement

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("COMPLEMENT"))

saved <- options(repr.plot.width = 5, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=c(3,6), targets.use=c(17,19), comparison=c(1,2), max.dataset=2, 
                        pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, line.on = TRUE,
                        line.size = 0.5, remove.isolate=F, dot.size.min = 5, dot.size.max = 10, font.size = 0, grid.on = TRUE,
                        color.grid = "grey75");
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

#limit to MIF

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("MIF"))

saved <- options(repr.plot.width = 5, repr.plot.height = 4.9) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=c(3,6), targets.use=c(17,19), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=F);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

#limit to MIF

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("MIF"))

saved <- options(repr.plot.width = 5, repr.plot.height = 4.9) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=c(3,6), targets.use=c(17), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=F);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

#limit to MIF

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("MIF"))

saved <- options(repr.plot.width = 5, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=c(3,6), targets.use=c(17), comparison=c(1,2), max.dataset=2, 
                        pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, line.on = TRUE,
                        line.size = 0.5, remove.isolate=F, dot.size.min = 5, dot.size.max = 14, font.size = 0, grid.on = TRUE,
                        color.grid = "grey75");
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(17:19), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(17:19), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 10, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(17,19), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(17,19), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 6, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=4, targets.use=c(17:19), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=4, targets.use=c(17:19), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 15) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:10), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:10), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

#Limit to TGFb

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("TGFb"))


saved <- options(repr.plot.width = 5.6, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:8), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=F);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:8), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=F);
gg1#+gg2;

options(saved) #restore old settings

#Limit to TGFb

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("TGFb"))


saved <- options(repr.plot.width = 5, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7,8), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, 
                        title.name="Increased signaling in SrfKO", angle.x=45, line.on = TRUE, line.size = 0.5, remove.isolate=F, 
                        dot.size.min = 5, dot.size.max = 14, font.size = 0, grid.on = TRUE,
                        color.grid = "grey75");
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:8), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=F);
gg1#+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 10, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:8), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=17, targets.use=c(7:8), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

#limit to IL1

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("IL1"))

saved <- options(repr.plot.width = 10, repr.plot.height = 5) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=c(17,19), targets.use=c(7), comparison=c(1,2), max.dataset=2, pairLR.use = pairLR.use, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=c(17,19), targets.use=c(7), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 6, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=18, targets.use=c(7:10), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=18, targets.use=c(7:10), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 15) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=19, targets.use=c(7:10), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=19, targets.use=c(7:10), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 6, repr.plot.height = 15) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=18, targets.use=c(3,4,6), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=18, targets.use=c(3,4,6), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(7:10), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=6, targets.use=c(7:10), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:10), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:10), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 10, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:8), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:8), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

#limit to EGF

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("EGF"))

saved <- options(repr.plot.width = 4.2, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
#gg1 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:8), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:8), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
#gg1+
gg2;

options(saved) #restore old settings

#limit to EGF

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("EGF"))

saved <- options(repr.plot.width = 5, repr.plot.height = 3) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 

gg2 <- netVisual_bubble(CCobjMerged, sources.use=3, targets.use=c(7:8), comparison=c(1,2), max.dataset=1, pairLR.use = pairLR.use, 
                        title.name="Decreased signaling in SrfKO", angle.x=45, line.on = TRUE, line.size = 0.5, remove.isolate=T, 
                        dot.size.min = 5, dot.size.max = 14, font.size = 0, grid.on = TRUE, color.grid = "grey75");
gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 6, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=4, targets.use=c(7:10), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
#gg2 <- netVisual_bubble(CCobjMerged, sources.use=4, targets.use=c(7:10), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1#+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(3,6), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(3,6), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

#limit to collagens
#see ?netVisual_bubble

saved <- options(repr.plot.width = 10, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

pairLR.use <- extractEnrichedLR(CCobjMerged, signaling = c("COLLAGEN"))

#CellChat can identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset. 
gg1 <- netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(3,6), comparison=c(1,2), max.dataset=2, title.name="Increased signaling in SrfKO", pairLR.use = pairLR.use, angle.x=45, remove.isolate=T);
gg2 <- netVisual_bubble(CCobjMerged, sources.use=7, targets.use=c(3,6), comparison=c(1,2), max.dataset=1, title.name="Decreased signaling in SrfKO", pairLR.use = pairLR.use, angle.x=45, remove.isolate=T);
gg1+gg2;

options(saved) #restore old settings

## alternative appoach: identify the up- and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis
CCobjMerged <- identifyOverExpressedGenes(CCobjMerged, group.dataset="datasets", pos.dataset="SrfKO", features.name="SrfKO", only.pos=FALSE, thresh.pc=0.1, thresh.fc=0.1, thresh.p=1);
net <- netMappingDEG(CCobjMerged, features.name="SrfKO");
net.up <- subsetCommunication(CCobjMerged, net=net, datasets="SrfKO", ligand.logFC=0.1, receptor.logFC=NULL);
net.down <- subsetCommunication(CCobjMerged, net=net, datasets="Ctrl", ligand.logFC=-0.1, receptor.logFC=-0.1);
gene.up <- extractGeneSubsetFromPair(net.up, CCobjMerged);
gene.down <- extractGeneSubsetFromPair(net.down, CCobjMerged);
pairLR.use.up <- net.up[, "interaction_name", drop=F];
pairLR.use.down <- net.down[, "interaction_name", drop=F];

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=7, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;

options(saved) #restore old settings

library(presto)
packageVersion("presto");

## alternative appoach: identify the up- and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis
CCobjMerged <- identifyOverExpressedGenes(CCobjMerged, group.dataset="datasets", pos.dataset="SrfKO", features.name="SrfKO", only.pos=FALSE, thresh.pc=0.1, thresh.fc=0.1, thresh.p=1);
net <- netMappingDEG(CCobjMerged, features.name="SrfKO");
net.up <- subsetCommunication(CCobjMerged, net=net, datasets="SrfKO", ligand.logFC=0.1, receptor.logFC=NULL);
net.down <- subsetCommunication(CCobjMerged, net=net, datasets="Ctrl", ligand.logFC=-0.1, receptor.logFC=-0.1);
gene.up <- extractGeneSubsetFromPair(net.up, CCobjMerged);
gene.down <- extractGeneSubsetFromPair(net.down, CCobjMerged);
pairLR.use.up <- net.up[, "interaction_name", drop=F];
pairLR.use.down <- net.down[, "interaction_name", drop=F];

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=7, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;

options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=7, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=7, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=6, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=6, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=3, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=3, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=4, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=4, targets.use=c(17:19), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=17, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=17, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=18, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=18, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 20) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=19, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=19, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=19, targets.use=c(3,4,6), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=19, targets.use=c(3,4,6), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=18, targets.use=c(3,4,6), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=18, targets.use=c(3,4,6), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=17, targets.use=c(3,4,6), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=17, targets.use=c(3,4,6), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=6, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=6, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=3, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=3, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

saved <- options(repr.plot.width = 12, repr.plot.height = 30) # Make the plots bigger from here on out, save old options

gg1 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.up, sources.use=4, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Up-regulated signaling in ", names(object.list)[2]));
gg2 <- netVisual_bubble(CCobjMerged, pairLR.use=pairLR.use.down, sources.use=4, targets.use=c(7:10), comparison=c(1,2), angle.x=90, remove.isolate=T, title.name=paste0("Down-regulated signaling in ", names(object.list)[2]));
gg1+gg2;
options(saved) #restore old settings

# Chord diagram - upregulated FibroF2 - Neutrophil

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 7, targets.use = c(19), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated Lum epi - macrophage

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 6, targets.use = c(17), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated Lum epi - neutrophil

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 6, targets.use = c(19), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated Lum epi - monocyte

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 6, targets.use = c(18), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated gland epi - macrophage

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 3, targets.use = c(17), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated gland epi - neutrophil

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 3, targets.use = c(19), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated gland epi - monocyte

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 3, targets.use = c(18), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated prolif epi - macrophage

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(17), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated prolif epi - neutrophil

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(19), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated prolif epi - monocyte

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(18), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated macrophage - fibroF2

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 17, targets.use = c(7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated neutrophil - fibroF2

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 19, targets.use = c(7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated monocyte - fibroF2

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 18, targets.use = c(7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - upregulated monocyte - lum ep

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 18, targets.use = c(6), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - lum ep - fibroF2

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 6, targets.use = c(7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - gland ep - fibroF2

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 3, targets.use = c(7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - gland ep - fibroF2

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(7), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

# Chord diagram - downregulated gland epi - FibroF2

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

netVisual_chord_gene(object.list[[1]], sources.use = 3, targets.use = c(7), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

options(saved) #restore old settings

## Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram.
## NOTE FROM SARA -- Most plots in part IV require specifying pathways and sometimes also the source & target cell types; you'll want to modify the signaling/sources.use/targets.use parameters according to your interests.

saved <- options(repr.plot.width = 12, repr.plot.height = 12) # Make the plots bigger from here on out, save old options

pathways.show <- c("CXCL"); 
weight.max <- getMaxWeight(object.list, slot.name=c("netP"), attribute=pathways.show);
par(mfrow=c(1,1), xpd=TRUE);
netVisual_aggregate(object.list[[1]], signaling=pathways.show, layout="circle", edge.weight.max=weight.max[1], edge.width.max=10, signaling.name=paste(pathways.show, names(object.list)[1]));

options(saved) #restore old settings

## Part IV: Visually compare cell-cell communication using Hierarchy plot, Circle plot or Chord diagram.
## NOTE FROM SARA -- Most plots in part IV require specifying pathways and sometimes also the source & target cell types; you'll want to modify the signaling/sources.use/targets.use parameters according to your interests.

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("MIF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("CCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("COLLAGEN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("COMPLEMENT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("CSF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("EGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("FGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("HGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("IFN-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("IL1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("IL17") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("MHC-I") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("MHC-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("WNT") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

pathways.show <- c("CXCL") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))

png("cellchat_outputs/merged.example_network_by_genotype.heatmap-view.png", h=1000, w=2000, res=200);
ht1 <- netVisual_heatmap(object.list[[1]], signaling=pathways.show, color.heatmap="Reds", title.name=paste(pathways.show, "signaling ", names(object.list)[1]));
ht2 <- netVisual_heatmap(object.list[[2]], signaling=pathways.show, color.heatmap="Reds", title.name=paste(pathways.show, "signaling ", names(object.list)[2]));
draw(ht1+ht2, ht_gap=unit(0.5, "cm"));
dev.off();

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

netVisual_heatmap(object.list[[1]], signaling=pathways.show, color.heatmap="Reds", title.name=paste(pathways.show, "signaling ", names(object.list)[1]));

options(saved) #restore old settings

# Chord diagram

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

options(saved) #restore old settings

pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow=c(1,1), xpd=TRUE);
netVisual_aggregate(object.list[[1]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list)[1]));


netVisual_aggregate(object.list[[2]], signaling=pathways.show, layout="chord", signaling.name=paste(pathways.show, names(object.list)[2]));

# (C) Identify signals contributing the most to outgoing or incoming signaling of certain cell groups

# In this heatmap, colobar represents the relative signaling strength of a signaling pathway across cell groups (NB: values are row-scaled).
# The top colored bar plot shows the total signaling strength of a cell group by summarizing all signaling pathways displayed in the heatmap.
# The right grey bar plot shows the total signaling strength of a signaling pathway by summarizing all cell groups displayed in the heatmap.

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(CCobjSrfKO, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(CCobjSrfKO, pattern = "incoming")

ht1
ht2

# Signaling role analysis on the cell-cell communication networks of interest in SrfKO only
netAnalysis_signalingRole_heatmap(CCobjSrfKO, signaling = c("COMPLEMENT", "IFN-II"))

# All the signaling pathways showing significant communications can be accessed by
CCobjSrfKO@netP$pathways

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CXCL") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CXCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CXCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CXCL") 
netAnalysis_contribution(CCobjCtrl, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CXCL <- extractEnrichedLR(CCobjCtrl, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjCtrl, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("COMPLEMENT") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COMPLEMENT <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COMPLEMENT[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COMPLEMENT <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COMPLEMENT[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COMPLEMENT <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COMPLEMENT[3,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COMPLEMENT <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COMPLEMENT[4,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("COLLAGEN") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[9,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.COLLAGEN <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.COLLAGEN[10,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("MMP") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.MMP <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MMP[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("IL1") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.IL1 <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IL1[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.IL1 <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IL1[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.IL1 <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.IL1[3,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("MIF") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.MIF <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MIF[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.MIF <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MIF[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CSF") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CSF <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CSF[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CSF <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CSF[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CCL") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[3,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[4,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[5,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[6,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[7,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[8,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[9,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[10,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[11,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[12,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[13,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[14,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[15,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[16,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[17,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[18,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[19,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[20,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[21,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[22,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[23,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CCL <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CCL[24,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("FN1") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.FN1 <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.FN1[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.FN1 <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.FN1[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("CX3C") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.CX3C <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CX3C[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
pathways.show <- c("TGFb") 
netAnalysis_contribution(CCobjSrfKO, signaling = pathways.show)
options(saved) # restore old settings

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# We can also visualize the cell-cell communication mediated by a single ligand-receptor pair. 
# We provide a function extractEnrichedLR to extract all the significant interactions (L-R pairs) and related signaling genes for a given signaling pathway.
pairLR.TGFb <- extractEnrichedLR(CCobjSrfKO, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[2,] # show one ligand-receptor pair (in this case #1)

# Circle plot
netVisual_individual(CCobjSrfKO, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(CCobjMerged) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
# Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net

#Export as csv
write.csv(df.net, "Merged_CellChat_seurat_df.csv")

# Extract the inferred cellular communication network as a data frame
df.net_SrfKO <- subsetCommunication(CCobjSrfKO) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
# Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net_SrfKO

#Export as csv
write.csv(df.net_SrfKO, "SrfKO_CellChat_seurat_df.csv")

# Extract the inferred cellular communication network as a data frame
df.net_CCobjCtrl <- subsetCommunication(CCobjCtrl) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
# Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net_CCobjCtrl

#Export as csv
write.csv(df.net_CCobjCtrl, "Ctrl_CellChat_seurat_df.csv")
