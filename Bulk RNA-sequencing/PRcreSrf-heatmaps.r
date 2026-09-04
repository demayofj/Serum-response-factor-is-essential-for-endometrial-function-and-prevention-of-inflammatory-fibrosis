# quick-start guide for plotting a DGE heatmap from a DGE table and a relative expression table via ComplexHeatmap
# full user-guide can be found here: https://jokergoo.github.io/ComplexHeatmap-reference/book/

library(ComplexHeatmap)
library(circlize)

# import normalized counts table (or whatever relative expression data)
# e.g. where the first column is gene name and the remaining columns are relative expression values for each sample in the experiment

norm.counts <- read.csv("250315_Heatmap.csv")
genes <- norm.counts$Gene # save gene names
norm.counts$Gene <- NULL
Annotations <- norm.counts$Annotation # save annotations
norm.counts$Annotation <- NULL
norm.counts <- as.matrix(norm.counts) # convert to numeric matrix
rownames(norm.counts) <- genes # make gene names the row names
rm(genes)
head(norm.counts)
tail(norm.counts)

# center and scale the rlog counts table (i.e. produce Z-score expression values across the samples, for each gene)
# note: t() will transpose the matrix
# so, here, we are transposing the matrix, center and scaling with scale(), then transposing back
# the reason for this is that scale() by default wants to scale each column vector (samples here), but we want to scale gene-wise

scaled.norm.counts <- t(scale(t(norm.counts)))

head(scaled.norm.counts)


# Export normalized counts file if desired

write.csv(scaled.norm.counts,file="scaled.norm.counts.csv") # keeps the rownames

set.seed(1) # set seed for reproducibility
Heatmap(scaled.norm.counts, 
        clustering_distance_rows = "euclidean", # this is the default method
        width = unit(10, "cm"), # width of the heatmap
        name="Normalized counts", # name of the legend i.e. colorbar
        col = colorRamp2(c(-2, -1, 0, 1, 2), 
                         c("#00007F", "blue", "white", "red", "#7F0000")), # this is a custom color scale
        show_row_names=TRUE,
        show_column_names=TRUE,
        column_title="Pgr cre/+ Srf f/f Uterus RNA-seq", # set your own title
        cluster_columns=TRUE)

#Cluster but do not reorder dendrograms
set.seed(1) # set seed for reproducibility
Heatmap(scaled.norm.counts, 
        clustering_distance_rows = "euclidean", # this is the default method
        width = unit(10, "cm"), # width of the heatmap
        name="Normalized count", # name of the legend i.e. colorbar
        col = colorRamp2(c(-2, -1, 0, 1, 2), 
                         c("#00007F", "blue", "white", "red", "#7F0000")), # this is a custom color scale
        show_row_names=TRUE,
        show_column_names=TRUE,
        column_title="Pgr cre/+ Srf f/f Uterus RNA-seq", # set your own title
        cluster_columns=TRUE,
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE)


#Cluster but do not reorder dendrograms AND split rows (and name them)
set.seed(1) # set seed for reproducibility
Heatmap(scaled.norm.counts, 
        clustering_distance_rows = "euclidean", # this is the default method
        width = unit(10, "cm"), # width of the heatmap
        name="Normalized count", # name of the legend i.e. colorbar
        col = colorRamp2(c(-2, -1, 0, 1, 2), 
                         c("#00007F", "blue", "white", "red", "#7F0000")), # this is a custom color scale
        show_row_names=TRUE,
        show_column_names=TRUE,
        column_title="Pgr cre/+ Srf f/f Uterus RNA-seq", # set your own title
        row_split = 2, # splits rows into the number indicated based on dendrogram
        row_title = c("Downregulated", "Upregulated"), # name the two sections we split it into
        cluster_columns=TRUE,
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE)


#Stop clustering rows and columns - keep original order 
set.seed(1) # set seed for reproducibility
Heatmap(scaled.norm.counts, 
        width = unit(10, "cm"), # width of the heatmap
        name="Normalized counts", # name of the legend i.e. colorbar
        col = colorRamp2(c(-2, -1, 0, 1, 2), 
                         c("#00007F", "blue", "white", "red", "#7F0000")), # this is a custom color scale
        show_row_names=TRUE,
        show_column_names=TRUE,
        column_title="Pgr cre/+ Srf f/f Uterus RNA-seq", # set your own title
        row_order = order(as.numeric(gsub("row", "", rownames(scaled.norm.counts)))), # orders rows by table order rather than clustering
        column_order = order(as.numeric(gsub("column", "", colnames(scaled.norm.counts))))) #orders columns by table order rather than clustering

# Adjust plot size
saved <- options(repr.plot.width = 8, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

#Stop clustering rows and columns
set.seed(1) # set seed for reproducibility
Heatmap(scaled.norm.counts, 
        width = unit(10, "cm"), # width of the heatmap
        name="Normalized counts", # name of the legend i.e. colorbar
        col = colorRamp2(c(-2, -1, 0, 1, 2), 
                         c("#00007F", "blue", "white", "red", "#7F0000")), # this is a custom color scale
        show_row_names=TRUE,
        show_column_names=TRUE,
        column_title="Pgr cre/+ Srf f/f Uterus RNA-seq", # set your own title
        row_order = order(as.numeric(gsub("row", "", rownames(scaled.norm.counts)))), # orders rows by table order rather than clustering
        column_order = order(as.numeric(gsub("column", "", colnames(scaled.norm.counts))))) #orders columns by table order rather than clustering

options(saved) # restore old settings

saved <- options(repr.plot.width = 8, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

# Change to default color scheme

# Stop clustering rows and columns
set.seed(1) # set seed for reproducibility
Heatmap(scaled.norm.counts, 
        width = unit(10, "cm"), # width of the heatmap
        name="Normalized counts", # name of the legend i.e. colorbar
        show_row_names=TRUE,
        show_column_names=TRUE,
        column_title="Pgr cre/+ Srf f/f Uterus RNA-seq", # set your own title
        row_order = order(as.numeric(gsub("row", "", rownames(scaled.norm.counts)))), # orders rows by table order rather than clustering
        column_order = order(as.numeric(gsub("column", "", colnames(scaled.norm.counts))))) #orders columns by table order rather than clustering

options(saved) # restore old settings

# FYI: to add back annotations
scaled.norm.counts.annotated <- cbind(scaled.norm.counts, Annotations)
head(scaled.norm.counts.annotated)

# Following: https://www.biostars.org/p/286187/#286507
# Load one more package
library(cluster)

# Load and view data
df <- read.csv("250315_Heatmap.csv")
head(df)

# center and scale the rlog counts table (i.e. produce Z-score expression values across the samples, for each gene)
# note: t() will transpose the matrix, and we are specifying that we are only taking the data starting in the 3rd column
# so, here, we are transposing the matrix, center and scaling with scale(), then transposing back
# the reason for this is that scale() by default wants to scale each column vector (samples here), but we want to scale gene-wise

heat <- t(scale(t(df[,3:ncol(df)])))
rownames(heat) <- df$Gene # assign row names as gene names
head(heat) # check our work

# Export normalized counts file if desired

write.csv(heat,file="heat_scaled.norm.counts.csv") # keeps the rownames

# Set annotations
  ColAnn <- data.frame(colnames(heat))
  colnames(ColAnn) <- c("Sample")
  ColAnn <- HeatmapAnnotation(df=ColAnn, which="col")

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c("Proinflammatory cytokines"="red","Interferon response"="coral", "Innate inflammatory response"="deeppink", "Myeloid immune cell markers"="darkred", "T cell markers"="darkmagenta", "B cell markers"="slateblue", "Estrogen signaling"="limegreen", "Progesterone signaling"="deepskyblue", "Cytoskeleton and contractility"="darkcyan", "Extracellular matrix"="darkseagreen", "Growth response and proliferation"="midnightblue", "Epithelial identity"="dimgrey"))
  RowAnn <- HeatmapAnnotation(df=RowAnn, col=colours, which="row")

myBreaks <- seq(-3, 3, length.out=100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)


# Set heatmap parameters and draw

saved <- options(repr.plot.width = 16, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

hmap <- Heatmap(heat,
        name="PRcre Srf RNA-seq",
        col=colorRamp2(myBreaks, myCol),
        heatmap_legend_param=list(color_bar="continuous", legend_direction="horizontal", legend_width=unit(5,"cm"), title_position="topcenter", title_gp=gpar(fontsize=15, fontface="bold")),

       #Split heatmap rows by annotation
        split=df$Annotation,

        #Row annotation configurations
        cluster_rows=TRUE,
        show_row_dend=TRUE,
        #row_title="Transcript", #overridden by 'split' it seems
        row_title_side="left",
        row_title_gp=gpar(fontsize=15, fontface="bold"),
        show_row_names=TRUE,
        row_names_side="left",
        row_title_rot=0,

        #Column annotation configuratiions
        cluster_columns=TRUE,
        show_column_dend=TRUE,
        column_title="Samples",
        column_title_side="top",
        column_title_gp=gpar(fontsize=15, fontface="bold"),
        column_title_rot=0,
        show_column_names=TRUE,

        #Dendrogram configurations: columns
        clustering_distance_columns=function(x) as.dist(1-cor(t(x))),
        clustering_method_columns="ward.D2",
        column_dend_height=unit(30,"mm"),

        #Dendrogram configurations: rows
        clustering_distance_rows="euclidean",
        clustering_method_rows="ward.D2",
        row_dend_width=unit(30,"mm"),

        #Annotations (row annotation must be added with 'draw' function, below)
        
        top_annotation=ColAnn)


draw(hmap + RowAnn, heatmap_legend_side="left", annotation_legend_side="right")

options(saved) # restore old settings

# Stop clustering samples

saved <- options(repr.plot.width = 16, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

hmap <- Heatmap(heat,
        name="PRcre Srf RNA-seq",
        col=colorRamp2(myBreaks, myCol),
        heatmap_legend_param=list(color_bar="continuous", legend_direction="horizontal", legend_width=unit(5,"cm"), title_position="topcenter", title_gp=gpar(fontsize=15, fontface="bold")),

       #Split heatmap rows by annotation
        split=df$Annotation,

        #Row annotation configurations
        cluster_rows=TRUE,
        show_row_dend=TRUE,
        #row_title="Transcript", #overridden by 'split' it seems
        row_title_side="left",
        row_title_gp=gpar(fontsize=15, fontface="bold"),
        show_row_names=TRUE,
        row_names_side="left",
        row_title_rot=0,

        #Column annotation configuratiions
        cluster_columns=FALSE,
        show_column_dend=FALSE,
        column_title="Samples",
        column_title_side="top",
        column_title_gp=gpar(fontsize=15, fontface="bold"),
        column_title_rot=0,
        show_column_names=TRUE,

        #Keep column order instead of clustering
        column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
                
        #Dendrogram configurations: rows
        clustering_distance_rows="euclidean",
        clustering_method_rows="ward.D2",
        row_dend_width=unit(30,"mm"),

        #Annotations (row annotation must be added with 'draw' function, below)
        
        top_annotation=ColAnn)


draw(hmap + RowAnn, heatmap_legend_side="left", annotation_legend_side="right")

options(saved) # restore old settings

# If you want to adjust the sample order but keep clustering turned off, change column order in the matrix
new_order <- c(1, 2, 4, 3, 8, 6, 7, 5)
heat <- heat[, new_order]
head(heat) # check our work

# Now plot again without sample annotation at the top

saved <- options(repr.plot.width = 16, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

#Set annotation parameters

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c("Proinflammatory cytokines"="red","Interferon response"="coral", "Innate inflammatory response"="deeppink", "Myeloid immune cell markers"="darkred", "T cell markers"="darkmagenta", "B cell markers"="slateblue", "Estrogen signaling"="limegreen", "Progesterone signaling"="deepskyblue", "Cytoskeleton and contractility"="darkcyan", "Extracellular matrix"="darkseagreen", "Growth response and proliferation"="midnightblue", "Epithelial identity"="dimgrey"))
  RowAnn <- HeatmapAnnotation(df=RowAnn, col=colours, which="row")

myBreaks <- seq(-3, 3, length.out=100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

#Set heatmap parameters
hmap <- Heatmap(heat,
        name="Normalized Counts",
        col=colorRamp2(myBreaks, myCol),
        heatmap_legend_param=list(color_bar="continuous", 
            legend_direction="horizontal", 
            legend_width=unit(5,"cm"), 
            title_position="topcenter", 
            title_gp=gpar(fontsize=14, fontface="bold"),
            labels_gp=gpar(fontsize=14)),

       #Split heatmap rows by annotation
        split=df$Annotation,

        #Row annotation configurations
        cluster_rows=TRUE,
        show_row_dend=TRUE,
        #row_title="Transcript", #overridden by 'split' it seems
        row_title_side="left",
        row_title_gp=gpar(fontsize=15, fontface="bold"),
        show_row_names=TRUE,
        row_names_side="left",
        row_title_rot=0,

        #Column annotation configuratiions
        cluster_columns=FALSE,
        show_column_dend=FALSE,
        column_title="", # blank the title and it won't show up
        column_title_side="top",
        column_title_gp=gpar(fontsize=15, fontface="bold"),
        column_title_rot=0,
        show_column_names=TRUE,

        #Keep column order instead of clustering
        column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
                
        #Dendrogram configurations: rows
        clustering_distance_rows="euclidean",
        clustering_method_rows="ward.D2",
        row_dend_width=unit(30,"mm"))

        #Annotations (row annotation must be added with 'draw' function, below)
        
        # top_annotation=ColAnn) comment out the top annotation


draw(hmap + RowAnn, heatmap_legend_side="left", annotation_legend_side="right")

options(saved) # restore old settings

# Move legend to the right side

saved <- options(repr.plot.width = 24, repr.plot.height = 24) # Make the plots bigger from here on out, save old options

set.seed(123)
#Set annotation parameters

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c("Proinflammatory cytokines"="red","Interferon response"="coral", "Innate inflammatory response"="deeppink", "Myeloid immune cell markers"="darkred", "T cell markers"="darkmagenta", "B cell markers"="slateblue", "Estrogen signaling"="limegreen", "Progesterone signaling"="deepskyblue", "Cytoskeleton and contractility"="darkcyan", "Extracellular matrix"="darkseagreen", "Growth response and proliferation"="midnightblue", "Epithelial identity"="dimgrey"))
  RowAnn <- HeatmapAnnotation(df=RowAnn, 
                              col=colours,
                              annotation_name_gp = gpar(fontsize = 16, fontface = "bold"),
                              which="row",
                              annotation_legend_param = list(
                                                title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                                labels_gp = gpar(fontsize = 16)))

myBreaks <- seq(-3, 3, length.out=100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

#Set heatmap parameters
hmap <- Heatmap(heat,
        name="Normalized Counts",
        col=colorRamp2(myBreaks, myCol),
        heatmap_width = unit(30,"cm"),
        heatmap_height = unit(50,"cm"),
        column_names_gp = gpar(fontsize = 16, fontface="bold"),
        heatmap_legend_param=list(color_bar="continuous", 
            legend_direction="horizontal", 
            legend_width=unit(7,"cm"), 
            title_position="topcenter", 
            title_gp=gpar(fontsize=16, fontface="bold"),
            labels_gp=gpar(fontsize=16)),

        #Keep column order instead of clustering
        column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
                
        #Split heatmap rows by annotation
        split=df$Annotation,

        #Row annotation configurations
        cluster_rows=TRUE,
        show_row_dend=TRUE,
        row_title_side="left",
        row_title_gp=gpar(fontsize=16, fontface="bold"),
        show_row_names=FALSE,
        row_names_side="left",
        row_title_rot=0,

        #Column annotation configuratiions
        cluster_columns=FALSE,
        show_column_dend=FALSE,
        column_title="",
        column_title_side="top",
        column_title_gp=gpar(fontsize=15, fontface="bold"),
        column_title_rot=0,
        show_column_names=TRUE,

        
                
        #Dendrogram configurations: rows
        clustering_distance_rows="euclidean",
        clustering_method_rows="ward.D2",
        row_dend_width=unit(10,"mm"))

        #Annotations (row annotation must be added with 'draw' function, below)
        


draw(hmap + RowAnn, heatmap_legend_side="right", annotation_legend_side="right")

options(saved) # restore old settings

# draw again without right side annotation (just left off the final draw function)

saved <- options(repr.plot.width = 24, repr.plot.height = 24) # Make the plots bigger from here on out, save old options

set.seed(123)
#Set annotation parameters

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c("Proinflammatory cytokines"="red","Interferon response"="coral", "Innate inflammatory response"="deeppink", "Myeloid immune cell markers"="darkred", "T cell markers"="darkmagenta", "B cell markers"="slateblue", "Estrogen signaling"="limegreen", "Progesterone signaling"="deepskyblue", "Cytoskeleton and contractility"="darkcyan", "Extracellular matrix"="darkseagreen", "Growth response and proliferation"="midnightblue", "Epithelial identity"="dimgrey"))
  RowAnn <- HeatmapAnnotation(df=RowAnn, 
                              col=colours,
                              annotation_name_gp = gpar(fontsize = 16, fontface = "bold"),
                              which="row",
                              annotation_legend_param = list(
                                                title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                                labels_gp = gpar(fontsize = 16)))

myBreaks <- seq(-3, 3, length.out=100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

#Set heatmap parameters
hmap <- Heatmap(heat,
        name="Normalized Counts",
        col=colorRamp2(myBreaks, myCol),
        heatmap_width = unit(30,"cm"),
        heatmap_height = unit(50,"cm"),
        column_names_gp = gpar(fontsize = 16, fontface="bold"),
        heatmap_legend_param=list(color_bar="continuous", 
            legend_direction="horizontal", 
            legend_width=unit(7,"cm"), 
            title_position="topcenter", 
            title_gp=gpar(fontsize=16, fontface="bold"),
            labels_gp=gpar(fontsize=16)),

        #Keep column order instead of clustering
        column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
                
        #Split heatmap rows by annotation
        split=df$Annotation,

        #Row annotation configurations
        cluster_rows=TRUE,
        show_row_dend=TRUE,
        row_title_side="left",
        row_title_gp=gpar(fontsize=16, fontface="bold"),
        show_row_names=FALSE,
        row_names_side="left",
        row_title_rot=0,

        #Column annotation configuratiions
        cluster_columns=FALSE,
        show_column_dend=FALSE,
        column_title="",
        column_title_side="top",
        column_title_gp=gpar(fontsize=15, fontface="bold"),
        column_title_rot=0,
        show_column_names=TRUE,

        
                
        #Dendrogram configurations: rows
        clustering_distance_rows="euclidean",
        clustering_method_rows="ward.D2",
        row_dend_width=unit(10,"mm"))

        #Annotations (row annotation must be added with 'draw' function, below)
        


draw(hmap, heatmap_legend_side="top")

options(saved) # restore old settings

# draw again without right side annotation, adding row titles

saved <- options(repr.plot.width = 24, repr.plot.height = 24) # Make the plots bigger from here on out, save old options

set.seed(123)
#Set annotation parameters

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c("Proinflammatory cytokines"="red","Interferon response"="coral", "Innate inflammatory response"="deeppink", "Myeloid immune cell markers"="darkred", "T cell markers"="darkmagenta", "B cell markers"="slateblue", "Estrogen signaling"="limegreen", "Progesterone signaling"="deepskyblue", "Cytoskeleton and contractility"="darkcyan", "Extracellular matrix"="darkseagreen", "Growth response and proliferation"="midnightblue", "Epithelial identity"="dimgrey"))
  RowAnn <- HeatmapAnnotation(df=RowAnn, 
                              col=colours,
                              annotation_name_gp = gpar(fontsize = 16, fontface = "bold"),
                              which="row",
                              annotation_legend_param = list(
                                                title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                                labels_gp = gpar(fontsize = 16)))

myBreaks <- seq(-3, 3, length.out=100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

#Set heatmap parameters
hmap <- Heatmap(heat,
        name="Normalized Counts",
        col=colorRamp2(myBreaks, myCol),
        heatmap_width = unit(30,"cm"),
        heatmap_height = unit(50,"cm"),
        column_names_gp = gpar(fontsize = 16, fontface="bold"),
        heatmap_legend_param=list(color_bar="continuous", 
            legend_direction="horizontal", 
            legend_width=unit(7,"cm"), 
            title_position="topcenter", 
            title_gp=gpar(fontsize=16, fontface="bold"),
            labels_gp=gpar(fontsize=16)),

        #Keep column order instead of clustering
        column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
                
        #Split heatmap rows by annotation
        split=df$Annotation,

        #Row annotation configurations
        cluster_rows=TRUE,
        show_row_dend=TRUE,
        row_title_side="left",
        row_title_gp=gpar(fontsize=16, fontface="bold"),
        show_row_names=TRUE,
        row_names_side="right",
        row_names_gp=gpar(fontsize=5.5),
        row_names_rot=0,
        row_title_rot=0,

        #Column annotation configuratiions
        cluster_columns=FALSE,
        show_column_dend=FALSE,
        column_title="",
        column_title_side="top",
        column_title_gp=gpar(fontsize=15, fontface="bold"),
        column_title_rot=0,
        show_column_names=TRUE,

        
                
        #Dendrogram configurations: rows
        clustering_distance_rows="euclidean",
        clustering_method_rows="ward.D2",
        row_dend_width=unit(10,"mm"))

        #Annotations (row annotation must be added with 'draw' function, below)
        


draw(hmap, heatmap_legend_side="top")

options(saved) # restore old settings

# draw again, now with right side annotation of specific genes
# See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation
# Also got some help from chatGPT to modify the code

saved <- options(repr.plot.width = 24, repr.plot.height = 32) # Make the plots bigger from here on out, save old options

set.seed(123)

# Generate heatmap colors
myBreaks <- seq(-3, 3, length.out = 100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

# Define rows to mark
rows_to_mark <- c("Il1f6", "Cxcl2", "Cxcl5", "Il1a", "Ccl5", "Ccl2", "Tnf",  # List genes to mark (from row names list); here proinflammatory cytokines
                  "Oas3", "Oas2", "Ifng", "Isg15", "Stat1", "Stat2", "Irf7", # interferon response
                  "Lcn2", "C3", "Nlrp3", "Zbp1", # innate inflammatory response
                  "Cxcr1", "Ly6g", "S100a9", "Cxcr2", "Ccr7", "Csf3r", "Fcgr4", # myeloid immune cell markers
                  "Ctla4", "Cd8a", "Cd4", # T cell markers
                  "Cd22", "Cd180", "Cd19", # B cell markers
                  "Sprr2d", "Ltf", "Clca3b", "Sprr2f", "Muc1", # estrogen signaling
                  "Zbtb16", "Scara5", "Fst", "Areg", "Ptch2", # progesterone signaling
                  "Srf", "Acta2", "Cnn1", "Tpm2", "Rhoj", "Des", "Vim", # cytoskeleton
                  "Mmp7", "Mmp8", "Col6a5", "Col12a1", "Col7a1", "Col17a1", "Vit", # extracellular matrix
                  "Egr1", "Fgf13", #growth response
                  "Wnt7a", "Epcam", "Prss28", "Spink1", "Prss29" # epithelial identity
                  )       
mark_indices <- match(rows_to_mark, rownames(heat))  # Get row indices

# Add row annotations using anno_mark
RowAnn <- rowAnnotation(
  mark = anno_mark(
    at = mark_indices,                 # Row indices to mark
    labels = rows_to_mark,             # Labels for the marked rows
    labels_gp = gpar(fontsize = 16, fontface = "italic"),   # Font size for labels
    link_gp = gpar(lwd = 1, col = "black") # Line properties
  )
)

# Set heatmap parameters
hmap <- Heatmap(
  heat,
  name = "Normalized Counts",
  col = colorRamp2(myBreaks, myCol),
  heatmap_width = unit(30, "cm"),
  heatmap_height = unit(50, "cm"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "horizontal",
    legend_width = unit(7, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  ),
  
    # Keep column order instead of clustering
  column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
  
    # Split heatmap rows by annotation
  split = df$Annotation,
  
    # Row annotation configurations
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 5.5),
  row_names_rot = 0,
  row_title_rot = 0,
  # Column annotation configurations
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = "",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  column_title_rot = 0,
  show_column_names = TRUE,
  # Dendrogram configurations: rows
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_width = unit(10, "mm")
)

# Draw heatmap with row annotations
draw(hmap + RowAnn, heatmap_legend_side = "top")

options(saved) # restore old settings


# Mark all genes, remove dendrogram

# draw again, now with right side annotation of specific genes
# See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation
# Also got some help from chatGPT to modify the code

saved <- options(repr.plot.width = 16, repr.plot.height = 42) # Make the plots bigger from here on out, save old options

set.seed(123)

# Generate heatmap colors
myBreaks <- seq(-3, 3, length.out = 100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

# Define rows to mark
rows_to_mark <- c("Il1f6", "Il1f9", "Cxcl2", "Cxcl5", "Il1a", "Il24", "Ccl20", "Ccl28", "Cxcl10", "Il1f8", "Ccl3", "Tnfsf15", "Cxcl11", "Ccl4", "Cxcl9", "Ccl5", "Cxcl1", "Cxcl17", "Cxcl3", "Tnf", "Ccl8", "Ccl12", "Ccl9", "Ccl2", "Cxcl16", "Ccl6", "Ccl22", "Tnfsf14", "Csf3", "Il7", "Oas3", "Ifit1bl1", "Ifng", "Oas2", "Ifne", "Ifit3b", "Oasl1", "Ifit3", "Ifi209", "Ifit1", "Ifi206", "Isg15", "Oas1g", "Ifi207", "Oas1a", "Oasl2", "Ifitm1", "Irf7", "Ifit1bl2", "Igtp", "Irf4", "Stat1", "Oas1h", "Ifit2", "Ifi204", "Ifi213", "Stat2", "Ifi44", "Oas1b", "Isg20", "Irf5", "Ifitm3", "Irf1", "Lcn2", "Wfdc21", "Wfdc18", "Zbp1", "C6", "C3", "Nlrp3", "Tlr13", "Tlr8", "Tlr7", "C2", "Nlrc5", "Tlr9", "C7", "Nlrc4", "Tlr2", "Tlr5", "Cd14", "Cxcr1", "Ly6g6g", "S100a9", "Cxcr2", "Ly6i", "Ly6c2", "Ccr7", "Cd300lf", "Ly6g", "Cd300lb", "Csf3r", "Fcgr4", "Fcgr1", "Lyz2", "Klra17", "Cd300ld", "Apoe", "C1qb", "C1qc", "C1qa", "Csf1r", "Ctla4", "Cd8a", "Cd3e", "Cxcr6", "Cd3g", "Cd3d", "Cd5", "Cd4", "Cd6", "Cd22", "Cd180", "Cd19", "Sprr2d", "Sprr2b", "Sprr2h", "Sprr2i", "Sprr2a2", "Ltf", "Clca3b", "Sprr2e", "Sprr2c-ps", "Clca1", "Sprr2f", "Wnt7b", "Sprr2g", "Wnt4", "Muc1", "Sgk1", "Klf4", "Hsd17b2", "Klf15", "Ptch1", "Ptch2", "Zbtb16", "Scara5", "Fst", "Areg", "Actg2", "Tagln", "Vim", "Mylk", "Rho", "Des", "Rhoj", "Myom1", "Cnn1", "Tpm1", "Tpm2", "Myl9", "Myh11", "Cnn2", "Acta2", "Srf", "Mmp12", "Mmp7", "Mmp13", "Mmp8", "Mmp10", "Col6a5", "Mmp27", "Col12a1", "Col6a6", "Mmp25", "Col7a1", "Col6a4", "Col8a2", "Col26a1", "Col8a1", "Mmp11", "Mmp19", "Timp3", "Col23a1", "Col15a1", "Col20a1", "Col17a1", "Mmp28", "Vit", "Fgf16", "Egfl6", "Egfl7", "Fgf13", "Egr1", "Egr3", "Wnt7a", "Epcam", "Foxj1", "Wnt9a", "Ttr", "Prss29", "Spink1", "Prss28", "Cdkn2b", "Cdkn1a", "Hgf", "Plcl1" 
                  )       
mark_indices <- match(rows_to_mark, rownames(heat))  # Get row indices

# Add row annotations using anno_mark
RowAnn <- rowAnnotation(
  mark = anno_mark(
    at = mark_indices,                 # Row indices to mark
    labels = rows_to_mark,             # Labels for the marked rows
    labels_gp = gpar(fontsize = 14, fontface = "italic"),   # Font size for labels
    link_gp = gpar(lwd = 1, col = "black") # Line properties
  )
)

# Set heatmap parameters
hmap <- Heatmap(
  heat,
  name = "Normalized Counts",
  col = colorRamp2(myBreaks, myCol),
  heatmap_width = unit(20, "cm"),
  heatmap_height = unit(100, "cm"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "horizontal",
    legend_width = unit(7, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  ),
  
    # Keep column order instead of clustering
  column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
  
    # Split heatmap rows by annotation
  split = df$Annotation,
  
    # Row annotation configurations
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 5.5),
  row_names_rot = 0,
  row_title_rot = 0,
  # Column annotation configurations
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = "",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  column_title_rot = 0,
  show_column_names = TRUE,
  # Dendrogram configurations: rows
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_width = unit(10, "mm")
)

# Draw heatmap with row annotations
draw(hmap + RowAnn, heatmap_legend_side = "top")

options(saved) # restore old settings


# draw again, now with right side annotation of specific genes
# See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation
# Also got some help from chatGPT to modify the code

saved <- options(repr.plot.width = 12, repr.plot.height = 21) # Make the plots bigger from here on out, save old options

set.seed(123)

# Generate heatmap colors
myBreaks <- seq(-3, 3, length.out = 100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

# Define rows to mark
rows_to_mark <- c("Il1f6", "Cxcl2", "Cxcl5", "Il1a", "Ccl5", "Ccl2", "Tnf",  # List genes to mark (from row names list); here proinflammatory cytokines
                  "Oas3", "Oas2", "Ifng", "Isg15", "Stat1", "Stat2", "Irf7", # interferon response
                  "Lcn2", "C3", "Nlrp3", "Zbp1", # innate inflammatory response
                  "Cxcr1", "Ly6g", "Cxcr2", "Ccr7", "Csf3r", "Fcgr4", # myeloid immune cell markers
                  "Ctla4", "Cd8a", "Cd4", # T cell markers
                  "Cd22", # B cell markers
                  "Sprr2d", "Ltf", "Clca3b", "Sprr2f", "Muc1", # estrogen signaling
                  "Zbtb16", "Scara5", "Fst", "Areg", # progesterone signaling
                  "Srf", "Acta2", "Cnn1", "Tpm2", "Des", "Vim", # cytoskeleton
                  "Mmp7", "Mmp8", "Col6a5", "Col12a1", "Col7a1", "Col17a1", "Vit", # extracellular matrix
                  "Egr1", "Fgf13", #growth response
                  "Epcam", "Prss28", "Spink1", "Prss29" # epithelial identity
                  )       
mark_indices <- match(rows_to_mark, rownames(heat))  # Get row indices

# Add row annotations using anno_mark
RowAnn <- rowAnnotation(
  mark = anno_mark(
    at = mark_indices,                 # Row indices to mark
    labels = rows_to_mark,             # Labels for the marked rows
    labels_gp = gpar(fontsize = 16, fontface = "italic"),   # Font size for labels
    link_gp = gpar(lwd = 1, col = "black") # Line properties
  )
)

# Set heatmap parameters
hmap <- Heatmap(
  heat,
  name = "Normalized Counts",
  col = colorRamp2(myBreaks, myCol),
  heatmap_width = unit(20, "cm"),
  heatmap_height = unit(50, "cm"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "horizontal",
    legend_width = unit(7, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  ),
  
    # Keep column order instead of clustering
  column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
  
    # Split heatmap rows by annotation
  split = df$Annotation,
  
    # Row annotation configurations
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 5.5),
  row_names_rot = 0,
  row_title_rot = 0,
  # Column annotation configurations
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = "",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  column_title_rot = 0,
  show_column_names = TRUE,
  # Dendrogram configurations: rows
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_width = unit(10, "mm")
)

# Draw heatmap with row annotations
draw(hmap + RowAnn, heatmap_legend_side = "top")

options(saved) # restore old settings


# draw again, now with right side annotation of specific genes
# See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation
# Also got some help from chatGPT to modify the code

saved <- options(repr.plot.width = 15, repr.plot.height = 21) # Make the plots bigger from here on out, save old options

set.seed(123)

# Generate heatmap colors
myBreaks <- seq(-3, 3, length.out = 100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c("Proinflammatory cytokines"="red","Interferon response"="coral", "Innate inflammatory response"="deeppink", "Myeloid immune cell markers"="darkred", "T cell markers"="darkmagenta", "B cell markers"="slateblue", "Estrogen signaling"="limegreen", "Progesterone signaling"="deepskyblue", "Cytoskeleton and contractility"="darkcyan", "Extracellular matrix"="darkseagreen", "Growth response and proliferation"="midnightblue", "Epithelial identity"="dimgrey"))
  RowAnn <- HeatmapAnnotation(df=RowAnn, 
                              col=colours,
                              annotation_name_gp = gpar(fontsize = 16, fontface = "bold"),
                              which="row",
                              annotation_legend_param = list(
                                                title_gp = gpar(fontsize = 16, fontface = "bold"), 
                                                labels_gp = gpar(fontsize = 16)))

# Set heatmap parameters
hmap <- Heatmap(
  heat,
  name = "Normalized Counts",
  col = colorRamp2(myBreaks, myCol),
  heatmap_width = unit(20, "cm"),
  heatmap_height = unit(50, "cm"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "horizontal",
    legend_width = unit(7, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  ),
  
    # Keep column order instead of clustering
  column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
  
    # Split heatmap rows by annotation
  split = df$Annotation,
  
    # Row annotation configurations
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 5.5),
  row_names_rot = 0,
  row_title_rot = 0,
  # Column annotation configurations
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = "",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  column_title_rot = 0,
  show_column_names = TRUE,
  # Dendrogram configurations: rows
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_width = unit(10, "mm")
)

# Draw heatmap with row annotations
draw(hmap + RowAnn, heatmap_legend_side = "right", annotation_legend_side="right")

options(saved) # restore old settings


# draw again, now with right side annotation of specific genes, move legend, modify gene selections
# See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation
# Also got some help from chatGPT to modify the code

saved <- options(repr.plot.width = 12, repr.plot.height = 16) # Make the plots bigger from here on out, save old options

set.seed(123)

# Generate heatmap colors
myBreaks <- seq(-3, 3, length.out = 100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

# Define rows to mark
rows_to_mark <- c("Il1f6", "Cxcl2", "Cxcl5", "Il1a", "Ccl5", "Ccl2", "Tnf", "Csf3", "Cxcl10",  # List genes to mark (from row names list); here proinflammatory cytokines
                  "Oas3", "Ifng", "Isg15", "Ifit3", "Oas1a", "Stat1", "Stat2", "Irf7", # interferon response
                  "Lcn2", "C3", "Nlrp3", "Zbp1", # innate inflammatory response
                  "Cxcr1", "Ly6g", "Cxcr2", "Ccr7", "Csf3r", "Fcgr4", # myeloid immune cell markers
                  "Ctla4", "Cd8a", "Cd4", # T cell markers
                  "Cd22", # B cell markers
                  "Sprr2d", "Ltf", "Clca3b", "Sprr2f", "Muc1", # estrogen signaling
                  "Scara5", "Fst", "Areg", # progesterone signaling
                  "Srf", "Acta2", "Cnn1", "Tpm2", "Des", "Vim", # cytoskeleton
                  "Mmp7", "Mmp8", "Col6a5", "Col12a1", "Col7a1", "Col17a1", "Vit", # extracellular matrix
                  "Egr1", "Fgf13", "Cdkn1a", #growth response
                  "Epcam", "Prss28", "Spink1" # epithelial identity
                  )       
mark_indices <- match(rows_to_mark, rownames(heat))  # Get row indices

# Add row annotations using anno_mark
RowAnn <- rowAnnotation(
  mark = anno_mark(
    at = mark_indices,                 # Row indices to mark
    labels = rows_to_mark,             # Labels for the marked rows
    labels_gp = gpar(fontsize = 16, fontface = "italic"),   # Font size for labels
    link_gp = gpar(lwd = 1, col = "black") # Line properties
  )
)

# Set heatmap parameters
hmap <- Heatmap(
  heat,
  name = "Normalized Counts",
  col = colorRamp2(myBreaks, myCol),
  heatmap_width = unit(20, "cm"),
  heatmap_height = unit(37, "cm"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(7, "cm"),
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  ),
  
    # Keep column order instead of clustering
  column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
  
    # Split heatmap rows by annotation
  split = df$Annotation,
  
    # Row annotation configurations
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 5.5),
  row_names_rot = 0,
  row_title_rot = 0,
  # Column annotation configurations
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = "",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  column_title_rot = 0,
  show_column_names = TRUE,
  # Dendrogram configurations: rows
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_width = unit(10, "mm")
)

# Draw heatmap with row annotations
draw(hmap + RowAnn, heatmap_legend_side = "right")

options(saved) # restore old settings


# Recategorize lymphoid immune markers

# Load and view data
df <- read.csv("250315_Heatmap.csv")
head(df)

# center and scale the rlog counts table (i.e. produce Z-score expression values across the samples, for each gene)
# note: t() will transpose the matrix, and we are specifying that we are only taking the data starting in the 3rd column
# so, here, we are transposing the matrix, center and scaling with scale(), then transposing back
# the reason for this is that scale() by default wants to scale each column vector (samples here), but we want to scale gene-wise

heat <- t(scale(t(df[,3:ncol(df)])))
rownames(heat) <- df$Gene # assign row names as gene names
head(heat) # check our work

# Export normalized counts file if desired

write.csv(heat,file="heat_scaled.norm.counts.csv") # keeps the rownames

# If you want to adjust the sample order but keep clustering turned off, change column order in the matrix
new_order <- c(1, 2, 4, 3, 8, 6, 7, 5)
heat <- heat[, new_order]
head(heat) # check our work

# draw again, now with right side annotation of specific genes, move legend, modify gene selections
# See https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#mark-annotation
# Also got some help from chatGPT to modify the code

saved <- options(repr.plot.width = 12, repr.plot.height = 16) # Make the plots bigger from here on out, save old options

set.seed(123)

# Generate heatmap colors
myBreaks <- seq(-3, 3, length.out = 100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)

# Define rows to mark
rows_to_mark <- c("Il1f6", "Cxcl2", "Cxcl5", "Il1a", "Ccl5", "Ccl2", "Tnf", "Csf3", "Cxcl10",  # List genes to mark (from row names list); here proinflammatory cytokines
                  "Oas3", "Ifng", "Isg15", "Ifit3", "Oas1a", "Stat1", "Stat2", "Irf7", # interferon response
                  "Lcn2", "C3", "Nlrp3", "Zbp1", # innate inflammatory response
                  "Cxcr1", "Ly6g", "Cxcr2", "Ccr7", "Csf3r", "Fcgr4", # myeloid immune cell markers
                  "Ctla4", "Cd8a", "Cd4", # T cell markers
                  "Cd22", # B cell markers
                  "Sprr2d", "Ltf", "Clca3b", "Sprr2f", "Muc1", # estrogen signaling
                  "Scara5", "Fst", "Areg", # progesterone signaling
                  "Srf", "Acta2", "Cnn1", "Tpm2", "Des", "Vim", # cytoskeleton
                  "Mmp7", "Mmp8", "Col6a5", "Col12a1", "Col7a1", "Col17a1", "Vit", # extracellular matrix
                  "Egr1", "Fgf13", "Cdkn1a", #growth response
                  "Epcam", "Prss28", "Spink1" # epithelial identity
                  )       
mark_indices <- match(rows_to_mark, rownames(heat))  # Get row indices

# Add row annotations using anno_mark
RowAnn <- rowAnnotation(
  mark = anno_mark(
    at = mark_indices,                 # Row indices to mark
    labels = rows_to_mark,             # Labels for the marked rows
    labels_gp = gpar(fontsize = 16, fontface = "italic"),   # Font size for labels
    link_gp = gpar(lwd = 1, col = "black") # Line properties
  )
)

# Set heatmap parameters
hmap <- Heatmap(
  heat,
  name = "Normalized Counts",
  col = colorRamp2(myBreaks, myCol),
  heatmap_width = unit(20, "cm"),
  heatmap_height = unit(37, "cm"),
  column_names_gp = gpar(fontsize = 16, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(7, "cm"),
    title_position = "leftcenter-rot",
    title_gp = gpar(fontsize = 16, fontface = "bold"),
    labels_gp = gpar(fontsize = 16)
  ),
  
    # Keep column order instead of clustering
  column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
  
    # Split heatmap rows by annotation
  split = df$Annotation,
  
    # Row annotation configurations
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 16, fontface = "bold"),
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 5.5),
  row_names_rot = 0,
  row_title_rot = 0,
  # Column annotation configurations
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  column_title = "",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 15, fontface = "bold"),
  column_title_rot = 0,
  show_column_names = TRUE,
  # Dendrogram configurations: rows
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "ward.D2",
  row_dend_width = unit(10, "mm")
)

# Draw heatmap with row annotations
draw(hmap + RowAnn, heatmap_legend_side = "right")

options(saved) # restore old settings

