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

# Read in the processed RDS file from Sara Grimm, which holds a seurat object with immune cells subsetted, integrated, and annotated (from 4 samples).
Adult_Srfflox_SrfKO_immune <- readRDS("/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-Immune.03feb2025.rds")
Adult_Srfflox_SrfKO_immune
head(Adult_Srfflox_SrfKO_immune)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="Sample", pt.size=0.5, reduction='umap', shuffle=TRUE);

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="genotype", pt.size=0.5, reduction='umap', shuffle=TRUE);

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="Phase", pt.size=0.5, reduction='umap', shuffle=TRUE)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="broad_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="broad_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="seurat_clusters", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       );

options(saved) # restore old settings

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE
       )+NoLegend();

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE
       )+NoLegend();

png(file = "Immune_PlotBySubtype.png",
    width = 16, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE
       )+NoLegend();

dev.off()

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=1, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype"
       )+NoLegend()

png(file = "Immune_PlotBySubtype_SplitByGenotype.png",
    width = 29, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=TRUE, split.by="genotype"
       )+NoLegend()

dev.off()

png(file = "Immune_PlotBySubtype_SplitByGenotype_nolabel.png",
    width = 29, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=0.5, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )+NoLegend()

dev.off()

png(file = "Immune_PlotBySubtype_SplitByGenotype_nolabel_largerpt.png",
    width = 29, height = 16, units = "cm", res = 600)

DimPlot(Adult_Srfflox_SrfKO_immune, group.by="ImmuneSubtype", pt.size=1, reduction='umap', shuffle=TRUE,
        label=FALSE, split.by="genotype"
       )+NoLegend()

dev.off()

# Use DittoSeq to plot fraction of each cluster by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_immune_tmp <- Adult_Srfflox_SrfKO_immune
Idents(Adult_Srfflox_SrfKO_immune_tmp) <- Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_immune_tmp)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
# Use dittoBarPlot to generate the plot with some customization

saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_immune_tmp, 
             var = "ImmuneSubtype", 
             group.by = "genotype",
             color.panel = my_color_palette, 
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings

# Use DittoSeq to plot fraction of each cluster by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_immune_tmp <- Adult_Srfflox_SrfKO_immune
Idents(Adult_Srfflox_SrfKO_immune_tmp) <- Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_immune_tmp)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))
# Use dittoBarPlot to generate the plot with some customization

saved <- options(repr.plot.width = 16, repr.plot.height = 10) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_immune_tmp, 
             var = "ImmuneSubtype", 
             group.by = "genotype",
             color.panel = my_color_palette,
             var.labels.reorder = c(6, 8, 3, 13, 12, 4, 7, 9, 11, 2, 10, 5, 1),
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 2/1)

options(saved) # restore old settings             


cluster_order

levels(Adult_Srfflox_SrfKO_immune_tmp)

my_color_palette

# Use DittoSeq to plot fraction of each cell type by genotype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

dittoBarPlot(Adult_Srfflox_SrfKO_immune_tmp, 
             var = "genotype", 
             group.by = "ImmuneSubtype", 
             main = "Cell Type Fractions") + theme(axis.text = element_text(size = 20), 
                    legend.text = element_text(size = 20), 
                    legend.key.height = unit(1, 'cm'), 
                    title = element_text(size = 20, hjust = 0.5), 
                    plot.title = element_text(hjust=0.5), 
                    aspect.ratio = 1/3)

options(saved) # restore old settings

# Count cells by type to compare between genotypes
cell_table <- table(Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype, Adult_Srfflox_SrfKO_immune_tmp$genotype)
cell_table

# Make a bar graph showing each immune subtype as a percentage of all cells by genotype (Calculated first in excel)
# help from ChIRP (GPT-4o)

# Create the dataset

ImmunePct_df <- data.frame(Immune_subtype = c("Macrophage", "Neutrophil", "cDC2_mregDC", "T_NKT", "ProliferatingNK", "gammadeltaT", "Monocyte", "NK", "Plasma", "cDC1", "pDC", "ILC2", "Bcell"),
  Percent_of_all_Ctrl_cells = c(2.665, 0.227, 0.602, 0.307, 1.222, 0.5, 0.051, 0.278, 0.381, 0.307, 0.119, 0.062, 0.028),
  Percent_of_all_SrfKO_cells = c(6.026, 4.824, 2.252, 1.749, 0.711, 1.377, 1.478, 1.044, 0.862, 0.528, 0.698, 0.113, 0.145))

# Print the dataframe
ImmunePct_df


# Reshape the data using tidyr
data_long <- ImmunePct_df %>%
  pivot_longer(cols = c(Percent_of_all_Ctrl_cells, Percent_of_all_SrfKO_cells),
               names_to = "Condition", values_to = "Percent")

data_long

# Adjust the Condition column for better readability
data_long$Condition <- gsub("Percent_of_all_", "", data_long$Condition)
data_long$Condition <- gsub("_cells", "", data_long$Condition)
data_long

# Create the bar plot
ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Use dittoseq colorscheme

ditto_colors <- dittoColors(get.names = FALSE) # This will return a vector of hex codes

# Create the bar plot
ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = ditto_colors) +
  theme_minimal() +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Create the bar plot with Seurat-like theme and dittoseq colorscheme
ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = ditto_colors)  # Custom colors similar to Seurat

# now match the colors to the seurat object

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_immune_tmp <- Adult_Srfflox_SrfKO_immune
Idents(Adult_Srfflox_SrfKO_immune_tmp) <- Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype
# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_immune_tmp)
# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

# Create the bar plot with Seurat-like theme and seurat object colorscheme
ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = my_color_palette)  # Custom colors similar to Seurat


#Split plot and color by genotype

# reorder data_long levels to match seurat object levels
data_long$Immune_subtype <- factor(data_long$Immune_subtype, levels=c("Macrophage", "Neutrophil", "cDC2_mregDC", "T_NKT", "ProliferatingNK", "gammadeltaT", "Monocyte", "NK", "Plasma", "cDC1", "pDC", "ILC2", "Bcell"))

ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Immune_subtype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = scales::hue_pal()(length(unique(data_long$Immune_subtype)))) +
  facet_wrap(~ Condition)

ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"))  # Customize colors for conditions

# now match the colors to the seurat object

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_immune_tmp <- Adult_Srfflox_SrfKO_immune
Idents(Adult_Srfflox_SrfKO_immune_tmp) <- Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype

# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_immune_tmp)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

# reorder data_long levels to match seurat object levels
data_long$Immune_subtype <- factor(data_long$Immune_subtype, levels=c("Macrophage", "Neutrophil", "cDC2_mregDC", "T_NKT", "ProliferatingNK", "gammadeltaT", "Monocyte", "NK", "Plasma", "cDC1", "pDC", "ILC2", "Bcell"))

#Plot
ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Immune_subtype)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) + 
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "top"
  ) +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = my_color_palette) +
  aes(group = Condition) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", linewidth = 0.5) #black outlines

# now match the colors to the seurat object and color outlines by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_immune_tmp <- Adult_Srfflox_SrfKO_immune
Idents(Adult_Srfflox_SrfKO_immune_tmp) <- Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype

# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_immune_tmp)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

# reorder data_long levels to match seurat object levels
data_long$Immune_subtype <- factor(data_long$Immune_subtype, levels=c("Macrophage", "Neutrophil", "cDC2_mregDC", "T_NKT", "ProliferatingNK", "gammadeltaT", "Monocyte", "NK", "Plasma", "cDC1", "pDC", "ILC2", "Bcell"))

# Define a custom color palette for conditions
condition_palette <- c('Ctrl' = '#949494', 'SrfKO' = '#000000')  # Gray for Ctrl and Black for SrfKO

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Immune_subtype, color = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), size = 1) + 
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.width = unit(2, "lines"),   # Increase the width of legend keys
    legend.key.height = unit(2, "lines"),  # Increase the height of legend keys
    legend.spacing.x = unit(0.4, "cm"),      # Increase horizontal spacing between legend items
    legend.spacing.y = unit(0.5, "cm"),       # Increase vertical spacing between legend items
    #aspect.ratio = 1/2
) +
  labs(
    title = "Percent of Cells by Immune Subtype",
    x = "Immune Subtype",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = my_color_palette) +
  scale_color_manual(values = condition_palette) +
  aes(group = Condition) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9))

options(saved) # restore old settings             


# flip
# now match the colors to the seurat object and color outlines by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_immune_tmp <- Adult_Srfflox_SrfKO_immune
Idents(Adult_Srfflox_SrfKO_immune_tmp) <- Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype

# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_immune_tmp)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

# reorder data_long levels to match seurat object levels
data_long$Immune_subtype <- factor(data_long$Immune_subtype, levels=c("Macrophage", "Neutrophil", "cDC2_mregDC", "T_NKT", "ProliferatingNK", "gammadeltaT", "Monocyte", "NK", "Plasma", "cDC1", "pDC", "ILC2", "Bcell"))

# Define a custom color palette for conditions
condition_palette <- c('Ctrl' = '#FFFFFF', 'SrfKO' = '#000000')  # White for Ctrl and Black for SrfKO

saved <- options(repr.plot.width = 10, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Immune_subtype, color = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), size = 0.8) + 
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(angle = 90, size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.width = unit(2, "lines"),   # Increase the width of legend keys
    legend.key.height = unit(2, "lines"),  # Increase the height of legend keys
    legend.spacing.x = unit(0.4, "cm"),      # Increase horizontal spacing between legend items
    legend.spacing.y = unit(0.5, "cm"),       # Increase vertical spacing between legend items
    #aspect.ratio = 1/2
) +
  labs(
    title = "",
    x = "",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = my_color_palette) +
  scale_color_manual(values = condition_palette) +
  aes(group = Condition) +
  geom_bar(stat = "identity", position = position_dodge(width = 1))

options(saved) # restore old settings             


# shorten
# now match the colors to the seurat object and color outlines by genotype

# first some prep work to extract the color pallette from the seurat object
require(scales)
Adult_Srfflox_SrfKO_immune_tmp <- Adult_Srfflox_SrfKO_immune
Idents(Adult_Srfflox_SrfKO_immune_tmp) <- Adult_Srfflox_SrfKO_immune_tmp$ImmuneSubtype

# Create vector with levels of object@ident
identities <- levels(Adult_Srfflox_SrfKO_immune_tmp)

# Create vector of default ggplot2 colors
my_color_palette <- hue_pal()(length(identities))

# reorder data_long levels to match seurat object levels
data_long$Immune_subtype <- factor(data_long$Immune_subtype, levels=c("Macrophage", "Neutrophil", "cDC2_mregDC", "T_NKT", "ProliferatingNK", "gammadeltaT", "Monocyte", "NK", "Plasma", "cDC1", "pDC", "ILC2", "Bcell"))

# Define a custom color palette for conditions
condition_palette <- c('Ctrl' = '#FFFFFF', 'SrfKO' = '#000000')  # White for Ctrl and Black for SrfKO

saved <- options(repr.plot.width = 10, repr.plot.height = 4.5) # Make the plots bigger from here on out, save old options

ggplot(data_long, aes(x = Immune_subtype, y = Percent, fill = Immune_subtype, color = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), size = 0.8) + 
  theme_classic() +  # Use a classic theme as a base
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(angle = 90, size = 10, color = "black"),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "right",
    legend.key.width = unit(2, "lines"),   # Increase the width of legend keys
    legend.key.height = unit(2, "lines"),  # Increase the height of legend keys
    legend.spacing.x = unit(0.4, "cm"),      # Increase horizontal spacing between legend items
    legend.spacing.y = unit(0.5, "cm"),       # Increase vertical spacing between legend items
    #aspect.ratio = 1/2
) +
  labs(
    title = "",
    x = "",
    y = "Percent of Cells"
  ) +
  scale_fill_manual(values = my_color_palette) +
  scale_color_manual(values = condition_palette) +
  aes(group = Condition) +
  geom_bar(stat = "identity", position = position_dodge(width = 1))

options(saved) # restore old settings             


head(Adult_Srfflox_SrfKO_immune_tmp)
tail(Adult_Srfflox_SrfKO_immune_tmp)

#Subset SrfKO only and Ctrl only
SrfKO_immune_tmp <- subset(Adult_Srfflox_SrfKO_immune_tmp, genotype == "SrfKO")

Ctrl_immune_tmp <- subset(Adult_Srfflox_SrfKO_immune_tmp, genotype == "Ctrl")

SrfKO_immune_tmp
Ctrl_immune_tmp


#plot immune markers

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('S100a9'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype",
            ) + theme(legend.position = "right")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('C1qa'), 
            pt.size=0.5,
            label=TRUE, 
            repel = TRUE,
            split.by="genotype",
            ) + theme(legend.position = "right")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Lyz2'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype",
            ) + theme(legend.position = "right")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Adgre1'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype",
            ) + theme(legend.position = "right")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Csf3r'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype",
            ) + theme(legend.position = "right")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('H2-Aa'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype",
            ) + theme(legend.position = "right")

options(saved) # restore old settings


#plot immune markers

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('S100a9'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('C1qa'), 
            pt.size=0.5,
            label=FALSE)
            
FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Lyz2'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Adgre1'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Csf3r'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Clec9a'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Il2rb'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Nkg7'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Cd3g'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Jchain'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ebf1'), 
            pt.size=0.5,
            label=FALSE)

options(saved) # restore old settings


#plot macrophage markers (see https://www.cellsignal.com/products/primary-antibodies/mouse-reactive-m1-vs-m2-macrophage-ihc-antibody-sampler-kit/97624)

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 8, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

# General

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Adgre1'), 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Cd68'), 
            pt.size=0.5,
            label=FALSE)

#M1

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Itgax'), #Cd11c
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Cd86'), 
            pt.size=0.5,
            label=FALSE)

#M2

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Mrc1'), #Cd206 
            pt.size=0.5,
            label=FALSE)

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Arg1'), 
            pt.size=0.5,
            label=FALSE)


options(saved) # restore old settings


# plot macrophage markers (see https://www.cellsignal.com/products/primary-antibodies/mouse-reactive-m1-vs-m2-macrophage-ihc-antibody-sampler-kit/97624)
# split by genotype and label

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

# General

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Adgre1'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Cd68'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")
#M1

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Itgax'), #Cd11c
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Cd86'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

#M2

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Mrc1'), #Cd206 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Arg1'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")


options(saved) # restore old settings


# plot more macrophage markers (see https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01084/full)
# split by genotype and label

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

#M1

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Acod1'), #Irg1
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Cd40'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ptgs2'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ccl2'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Csf2'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Irf7'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ifit2'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Irf9'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ifi35'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ifnar2'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Isg20'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ifih1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Il12a'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Il12b'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Jak2'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Tlr4'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Myd88'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Nfkb1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Rela'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings


# plot more macrophage markers (see https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01084/full)
# split by genotype and label

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

#M2

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Lpxn'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Dhrs3'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Mical1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Dnmt3a'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Jun'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Gab1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('P2ry1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Fbxo32'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Cebpa'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Gadd45g'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('N4bp2l1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Tns1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Lpin1'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Mnt'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Rgs2'),
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings


# plot monocyte/neutrophil markers - better done as surface labeling than gene expression
# split by genotype and label

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

# General

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ly6g'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Ly6c1'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings


# plot IL1 genes
# split by genotype and label

Idents(Adult_Srfflox_SrfKO_immune) <- Adult_Srfflox_SrfKO_immune$ImmuneSubtype

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options

# General

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Il1a'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

FeaturePlot(Adult_Srfflox_SrfKO_immune, features = c('Il1b'), 
            pt.size=0.5,
            label=TRUE,
            repel = TRUE,
            split.by="genotype")

options(saved) # restore old settings

