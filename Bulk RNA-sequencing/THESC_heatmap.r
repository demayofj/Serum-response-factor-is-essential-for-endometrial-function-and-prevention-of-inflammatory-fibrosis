# quick-start guide for plotting a DGE heatmap from a DGE table and a relative expression table via ComplexHeatmap
# full user-guide can be found here: https://jokergoo.github.io/ComplexHeatmap-reference/book/

library(ComplexHeatmap)
library(circlize)

# Load and view THESC normalized count data
df <- read.csv("250329_VehHeatmap_THESC_final.csv")
head(df)


# center and scale the rlog counts table (i.e. produce Z-score expression values across the samples, for each gene)
# note: t() will transpose the matrix, and we are specifying that we are only taking the data starting in the 3rd column
# so, here, we are transposing the matrix, center and scaling with scale(), then transposing back
# the reason for this is that scale() by default wants to scale each column vector (samples here), but we want to scale gene-wise

heat <- t(scale(t(df[,3:ncol(df)])))
rownames(heat) <- df$Gene # assign row names as gene names
head(heat) # check our work

# Export scaled normalized counts file if desired

write.csv(heat,file="250329_heat_scaled.norm.counts_THESC_final.csv") # keeps the rownames


# Set annotations
  ColAnn <- data.frame(colnames(heat))
  colnames(ColAnn) <- c("Sample")
  ColAnn <- HeatmapAnnotation(df=ColAnn, which="col")

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c("Extracellular matrix"="yellow", "Growth response and proliferation"="slateblue", "Cytoskeleton and contractility"="deepskyblue"))
  RowAnn <- HeatmapAnnotation(df=RowAnn, col=colours, which="row")

myBreaks <- seq(-3, 3, length.out=100)
myCol <- colorRampPalette(c("blue", "white", "red"))(100)


# draw without right side annotation, adding all row titles instead

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

set.seed(123)
#Set annotation parameters

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c( "Extracellular matrix"="yellow", "Growth response and proliferation"="slateblue", "Cytoskeleton and contractility"="deepskyblue"))
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
        heatmap_width = unit(18,"cm"),
        heatmap_height = unit(18,"cm"),
        column_names_gp = gpar(fontsize = 16, fontface="bold"),
        heatmap_legend_param = list(
            color_bar = "continuous",
            legend_direction = "vertical",
            legend_width = unit(7, "cm"),
            title_position = "leftcenter-rot",
            title_gp = gpar(fontsize = 16, fontface = "bold"),
            labels_gp = gpar(fontsize = 16)),

        #Keep row order instead of clustering
        #row_order = order(as.numeric(gsub("row", "", rownames(heat)))),
        #Keep column order instead of clustering
        #column_order = order(as.numeric(gsub("column", "", colnames(heat)))),
                
        #Split heatmap rows by annotation
        split=df$Annotation,

        #Row annotation configurations
        cluster_rows=FALSE,
        show_row_dend=FALSE,
        row_title_side="left",
        row_title_gp=gpar(fontsize=16, fontface="bold"),
        show_row_names=TRUE,
        row_names_side="right",
        row_names_gp=gpar(fontsize=12),
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
        #clustering_distance_rows="euclidean",
        #clustering_method_rows="ward.D2",
        #row_dend_width=unit(10,"mm")
        )

        #Annotations (row annotation must be added with 'draw' function, below)
        


draw(hmap, heatmap_legend_side="right")

options(saved) # restore old settings

#another version with rotated legend
# draw without right side annotation, adding all row titles instead

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

set.seed(123)
#Set annotation parameters

RowAnn <- data.frame(df$Annotation)
  colnames(RowAnn) <- c("Annotation")
  colours <- list("Annotation"=c( "Extracellular matrix"="yellow", "Growth response and proliferation"="slateblue", "Cytoskeleton and contractility"="deepskyblue"))
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
        heatmap_width = unit(18,"cm"),
        heatmap_height = unit(18,"cm"),
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
        cluster_rows=FALSE,
        show_row_dend=FALSE,
        row_title_side="left",
        row_title_gp=gpar(fontsize=16, fontface="bold"),
        show_row_names=TRUE,
        row_names_side="right",
        row_names_gp=gpar(fontsize=12),
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
        #clustering_distance_rows="euclidean",
        #clustering_method_rows="ward.D2",
        #row_dend_width=unit(10,"mm")
        )

        #Annotations (row annotation must be added with 'draw' function, below)
        


draw(hmap, heatmap_legend_side="top")

options(saved) # restore old settings
