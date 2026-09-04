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
library(stringr)
library(ComplexHeatmap)

# Record active package versions
packageVersion("Seurat");
packageVersion("ggplot2");        
packageVersion("dplyr");
packageVersion("viridis");
packageVersion("dittoSeq");
packageVersion("scales");
packageVersion("scCustomize");
packageVersion("tidyr");
packageVersion("stringr");
packageVersion("ComplexHeatmap");

# Read in the processed RDS file from Sara Grimm, which is the human endometrial cell atlas non-hormone treated (HECA-NH) seurat object.
# See prep code: "step0-make_NonHormoneDonor_HECA_object.Rtxt"
hecaNH <- readRDS("/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/hecaNH_integration.31mar2025.rds")

hecaNH
head(hecaNH)

levels(hecaNH$Binary_Stage)

table(hecaNH$binary_pathology,hecaNH$Binary_Stage);

table(hecaNH$binary_pathology,hecaNH$celltypeCollapsed);

table(hecaNH$lineage,hecaNH$Binary_Stage);

table(hecaNH$Binary_Stage,hecaNH$celltypeCollapsed);

# Plot genes of interest

saved <- options(repr.plot.width = 10, repr.plot.height = 10.5) # Make the plots bigger from here on out, save old options
FeaturePlot(hecaNH,
            features = c('SRF'),
            label = TRUE, order = FALSE, min.cutoff = 0, max.cutoff = 1.5,
            pt.size = 2)  + theme(legend.position = "right")
options(saved) # restore old settings


# Plot genes of interest

saved <- options(repr.plot.width = 10, repr.plot.height = 10.5) # Make the plots bigger from here on out, save old options
FeaturePlot(hecaNH,
            features = c('SRF'),
            label = FALSE, order = FALSE, min.cutoff = 0, max.cutoff = 1.5,
            pt.size = 2)  + theme(legend.position = "right")
options(saved) # restore old settings


# Plot genes of interest

saved <- options(repr.plot.width = 20, repr.plot.height = 10.5) # Make the plots bigger from here on out, save old options
FeaturePlot(hecaNH,
            features = c('SRF'),
            split.by = "binary_pathology",
            label = TRUE, order = FALSE, min.cutoff = 0, max.cutoff = 1.5,
            pt.size = 2)  + theme(legend.position = "right")
options(saved) # restore old settings


# Plot genes of interest

saved <- options(repr.plot.width = 20, repr.plot.height = 10.5) # Make the plots bigger from here on out, save old options
FeaturePlot(hecaNH,
            features = c('SRF'),
            split.by = "binary_pathology",
            label = FALSE, order = FALSE, min.cutoff = 0, max.cutoff = 1.5,
            pt.size = 2)  + theme(legend.position = "right")
options(saved) # restore old settings


# SRF violin Plot by cell type

Idents(hecaNH) <- hecaNH$celltypeCollapsed

saved <- options(repr.plot.width = 16, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
VlnPlot(hecaNH, features = c("SRF"), 
        group.by = "celltypeCollapsed",
        pt.size = 0.1)

options(saved) # restore old settings

# SRF violin Plot by cell type/pathology

Idents(hecaNH) <- hecaNH$celltypeCollapsed

saved <- options(repr.plot.width = 16, repr.plot.height = 8) # Make the plots bigger from here on out, save old options
VlnPlot(hecaNH, features = c("SRF"), 
        group.by = "celltypeCollapsed",
        pt.size = 0.1,
        cols = c("dimgrey", "firebrick2"),
        split.plot = TRUE,
        split.by="binary_pathology",
)

options(saved) # restore old settings

# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Create a new column for combined cell type and pathology
data$combined <- interaction(data$celltypeCollapsed, data$binary_pathology, sep = "_")

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_point(aes(color = binary_pathology), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 0.3) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
p


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "Other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Convert celltypeCollapsed back to factor
data$celltypeCollapsed <- factor(data$celltypeCollapsed)

# Reorder the levels of celltypeCollapsed to move "Endo" to the end
levels_order <- c(setdiff(levels(data$celltypeCollapsed), "Endo"), "Endo")
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = levels_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_jitter(aes(color = binary_pathology), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.1) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings

# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_jitter(aes(color = binary_pathology), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.3) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_jitter(aes(color = binary_pathology), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 16) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create a custom labeling function for y-axis
custom_y_labels <- function(breaks) {
  labels <- ifelse(breaks %% 1 == 0, as.character(breaks), as.character(breaks))
  return(labels)
}

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  scale_y_continuous(
    breaks = seq(floor(min(data$SRF)), ceiling(max(data$SRF)), by = 1),  # Define major breaks
    minor_breaks = seq(floor(min(data$SRF)), ceiling(max(data$SRF)), by = 0.5),  # Define minor breaks
    labels = custom_y_labels
  ) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 16) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create a custom labeling function for y-axis
custom_y_labels <- function(breaks) {
  labels <- ifelse(breaks %% 1 == 0, as.character(breaks), as.character(breaks))
  return(labels)
}

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  scale_y_continuous(
    breaks = seq(floor(min(data$SRF)), ceiling(max(data$SRF)), by = 1),  # Define major breaks
    minor_breaks = seq(floor(min(data$SRF)), ceiling(max(data$SRF)), by = 0.5),  # Define minor breaks
    labels = custom_y_labels
  ) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 3) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.8),
    colour = "black"
  ) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_jitter(aes(color = binary_pathology), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.3) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.8),
    colour = "black"
  ) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Determine y-axis limits based on the data, with some padding
y_max <- max(data$SRF, na.rm = TRUE)
y_min <- min(data$SRF, na.rm = TRUE)
padding <- 0.05 * (y_max - y_min)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.8),
    colour = "black"
  ) +
  geom_jitter(
    aes(color = binary_pathology),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 0.5,
    alpha = 0.6
  ) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(y_min - padding, y_max + padding))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "binary_pathology"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Determine y-axis limits based on the data, with some padding
y_max <- max(data$SRF, na.rm = TRUE)
y_min <- min(data$SRF, na.rm = TRUE)
padding <- 0.05 * (y_max - y_min)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = binary_pathology)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.8),
    colour = "black"
  ) +
  geom_jitter(
    aes(color = binary_pathology),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 0.5,
    alpha = 0.6
  ) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(y_min - padding, y_min + padding + 0.5))  # Adjust y-axis limits to zoom in

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "Binary_Stage"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Convert Binary_Stage to character for easier manipulation
data$Binary_Stage <- as.character(data$Binary_Stage)

# Ensure that desired groups are included in desired order
desired_order <- c("Proliferative", "Secretory", "Menstrual")
data <- data[data$Binary_Stage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$Binary_Stage <- factor(data$Binary_Stage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = Binary_Stage)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_jitter(aes(color = Binary_Stage), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.3) +
  scale_fill_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  scale_color_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "Binary_Stage"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Convert Binary_Stage to character for easier manipulation
data$Binary_Stage <- as.character(data$Binary_Stage)

# Ensure that desired groups are included in desired order
desired_order <- c("Proliferative", "Secretory", "Menstrual")
data <- data[data$Binary_Stage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$Binary_Stage <- factor(data$Binary_Stage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = Binary_Stage)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_jitter(aes(color = Binary_Stage), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  scale_color_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "Binary_Stage"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Convert Binary_Stage to character for easier manipulation
data$Binary_Stage <- as.character(data$Binary_Stage)

# Ensure that desired groups are included in desired order
desired_order <- c("Proliferative", "Secretory", "Menstrual")
data <- data[data$Binary_Stage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$Binary_Stage <- factor(data$Binary_Stage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = Binary_Stage)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  scale_fill_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  scale_color_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "Binary_Stage"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Convert Binary_Stage to character for easier manipulation
data$Binary_Stage <- as.character(data$Binary_Stage)

# Ensure that desired groups are included in desired order
desired_order <- c("Proliferative", "Secretory", "Menstrual")
data <- data[data$Binary_Stage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$Binary_Stage <- factor(data$Binary_Stage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = Binary_Stage)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  scale_color_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "Binary_Stage"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Convert Binary_Stage to character for easier manipulation
data$Binary_Stage <- as.character(data$Binary_Stage)

# Ensure that desired groups are included in desired order
desired_order <- c("Proliferative", "Secretory", "Menstrual")
data <- data[data$Binary_Stage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$Binary_Stage <- factor(data$Binary_Stage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = Binary_Stage)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  scale_color_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 12) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "Binary_Stage"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Convert Binary_Stage to character for easier manipulation
data$Binary_Stage <- as.character(data$Binary_Stage)

# Ensure that desired groups are included in desired order
desired_order <- c("Proliferative", "Secretory", "Menstrual")
data <- data[data$Binary_Stage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$Binary_Stage <- factor(data$Binary_Stage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = Binary_Stage)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, position = position_dodge(width = 0.8), colour = "black") +
  scale_fill_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  scale_color_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 3) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$celltypeCollapsed

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "celltypeCollapsed", "Binary_Stage"))

# Convert celltypeCollapsed to character for easier manipulation
data$celltypeCollapsed <- as.character(data$celltypeCollapsed)

# Filter out the "SMC" and "EndoLymph" groups because they have fewer that 50 cells in the endometriosis datasets, also "other"
data <- data[!data$celltypeCollapsed %in% c("SMC", "EndoLymph", "other"), ]

# Combine the "FibF3" and "FibF2" groups into a single group
data$celltypeCollapsed[data$celltypeCollapsed %in% c("FibF3", "FibF2")] <- "FibF"

# Ensure that only desired groups are included
desired_order <- c("EpiLum", "EpiGlnd", "EpiCyc", "FibF", "FibCyc", "Periv", "Endo", "ImmMyel", "ImmLymph")
data <- data[data$celltypeCollapsed %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$celltypeCollapsed <- factor(data$celltypeCollapsed, levels = desired_order)

# Convert Binary_Stage to character for easier manipulation
data$Binary_Stage <- as.character(data$Binary_Stage)

# Ensure that desired groups are included in desired order
desired_order <- c("Proliferative", "Secretory", "Menstrual")
data <- data[data$Binary_Stage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$Binary_Stage <- factor(data$Binary_Stage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = celltypeCollapsed, y = SRF, fill = Binary_Stage)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.8),
    colour = "black"
  ) +
  scale_fill_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  scale_color_manual(values = c("dodgerblue3", "chartreuse3", "darkorange2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 8, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# SRF violin Plot by cell type/pathology

Idents(hecaNH) <- hecaNH$celltypeCollapsed

saved <- options(repr.plot.width = 16, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
dittoPlot(hecaNH, "SRF", 
        group.by = "celltypeCollapsed",
        plots = c("jitter", "vlnplot"),
        jitter.color = c("dimgrey"), jitter.size = 0.1,
        split.by="binary_pathology",
        split.nrow = 2
)

options(saved) # restore old settings

# SRF dotplot by cell type

genes <- c("SRF")

Idents(hecaNH) <- hecaNH$celltypeCollapsed

saved <- options(repr.plot.width = 4, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

DotPlot(hecaNH, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');
       

options(saved) # restore old settings

# SRF violin Plot by cell lineage

Idents(hecaNH) <- hecaNH$lineage

saved <- options(repr.plot.width = 16, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
VlnPlot(hecaNH, features = c("SRF"), 
        group.by = "lineage",
        pt.size = 0.1)

options(saved) # restore old settings

# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$lineage

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "lineage", "binary_pathology"))

# Convert lineage to character for easier manipulation
data$lineage <- as.character(data$lineage)

# Ensure that only desired groups are included
desired_order <- c("Epithelial", "Mesenchymal", "Endothelial", "Immune")
data <- data[data$lineage %in% desired_order, ]

# Convert celltypeCollapsed back to factor with the desired order
data$lineage <- factor(data$lineage, levels = desired_order)

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = lineage, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_jitter(aes(color = binary_pathology), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.3) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Cell Type", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 6, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# SRF dotplot by cell lineage

genes <- c("SRF")

Idents(hecaNH) <- hecaNH$lineage

saved <- options(repr.plot.width = 4, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

DotPlot(hecaNH, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');
       

options(saved) # restore old settings

# SRF violin Plot by cycle stage

Idents(hecaNH) <- hecaNH$Binary_Stage

saved <- options(repr.plot.width = 16, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
VlnPlot(hecaNH, features = c("SRF"), 
        group.by = "Binary_Stage",
        pt.size = 0.1)

options(saved) # restore old settings

# Set the identities in the Seurat object
Idents(hecaNH) <- hecaNH$Binary_Stage

# Extract the data used to create the plot
data <- FetchData(hecaNH, vars = c("SRF", "Binary_Stage", "binary_pathology"))

# Create the base plot with ggplot2 to control the point colors
p <- ggplot(data, aes(x = Binary_Stage, y = SRF, fill = binary_pathology)) +
  geom_violin(position = position_dodge(width = 0.8), scale = "width", trim = FALSE) +
  geom_jitter(aes(color = binary_pathology), position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8), size = 0.3) +
  scale_fill_manual(values = c("dimgrey", "firebrick2")) +
  scale_color_manual(values = c("dimgrey", "firebrick2")) +
  theme_minimal() +
  labs(x = "Binary_Stage", y = "SRF Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
saved <- options(repr.plot.width = 4, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
p
options(saved) # restore old settings


# SRF dotplot by cycle stage

genes <- c("SRF")

Idents(hecaNH) <- hecaNH$Binary_Stage

saved <- options(repr.plot.width = 4, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

DotPlot(hecaNH, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');
       

options(saved) # restore old settings

# SRF violin Plot by pathology status

Idents(hecaNH) <- hecaNH$binary_pathology

saved <- options(repr.plot.width = 16, repr.plot.height = 4) # Make the plots bigger from here on out, save old options
VlnPlot(hecaNH, features = c("SRF"), 
        group.by = "binary_pathology",
        pt.size = 0.1)

options(saved) # restore old settings

# SRF dotplot by pathology status

genes <- c("SRF")

Idents(hecaNH) <- hecaNH$binary_pathology

saved <- options(repr.plot.width = 4, repr.plot.height = 4) # Make the plots bigger from here on out, save old options

DotPlot(hecaNH, features=genes) + RotatedAxis() + scale_color_viridis(option="D") + theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + ylab('') + xlab('');
       

options(saved) # restore old settings

head(hecaNH)

levels(hecaNH$celltypeCollapsed)

# Let's make a cell type and phase-specific dotplot follwing instructions from https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/

# Customized multi-group dotplot
# We need to get the percentage of positive cells and average expression by group.
# For a single gene, put the groups into multiple rows, and each column is a cell type.

# group1 is the cell type/cluster annotation 
# group2 is any condition you want to further group, in this case, the phase

#This first part is just copy and paste from the link to prepare the function; no customization

GetMatrixFromSeuratByGroupSingle<- function(obj, feature, group1, group2){
  if (length(feature) != 1){
          stop("please only provide only one gene name")
  }
  exp_mat<- obj@assays$RNA$data[feature, ,drop=FALSE]
  count_mat<- obj@assays$RNA$counts[feature,,drop=FALSE ]
  
  meta<- obj@meta.data %>%
  tibble::rownames_to_column(var = "cell")
        
  # get the average expression matrix
  exp_df<- as.matrix(exp_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta) %>%
    group_by(gene,{{group1}}, {{group2}}) %>%
    summarise(average_expression = mean(expression)) %>%
    tidyr::pivot_wider(names_from = {{group1}}, 
                       values_from= average_expression) 
  
  exp_mat<- exp_df[, -c(1,2)] %>% as.matrix()
  rownames(exp_mat)<- exp_df %>% pull({{group2}})
  
  # get the percentage positive cell matrix
  count_df<- as.matrix(count_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    left_join(meta) %>%
    group_by(gene, {{group1}}, {{group2}}) %>%
    summarise(percentage = mean(count >0)) %>%
    tidyr::pivot_wider(names_from = {{group1}}, 
                       values_from= percentage) 

  percent_mat<- count_df[, -c(1,2)] %>% as.matrix()
  rownames(percent_mat)<- count_df %>% pull({{group2}})
  
  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }
  
  if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  
  if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}



# Let’s get the matrices for one gene

# This is where we enter our object name, feature of interest, and condition names from the object metadata

mat<- GetMatrixFromSeuratByGroupSingle(obj= hecaNH, 
                                 feature = "SRF", 
                                 lineage,
                                 Binary_Stage)

# take a look at the matrices

# 1. the average expression for each cell type per condition
mat$exp_mat

# 2. the percentage of cells positive for the gene of interest for each cell type per condition
mat$percent_mat

# Now, Let’s visualize it using ComplexHeatmap
# Always explore the data range before you decide how to map the data values to colors.

quantile(mat$exp_mat, c(0.1, 0.5, 0.8, 0.9))

# In this case, 0 will be mapped to #FDE725FF, 0.126 will be mapped to #238A8DFF and 0.17 will be mapped to #440154FF. The value in-between will be linearlly interpolated to get corresponding colors

col_fun<-  circlize::colorRamp2(c(0, 0.126, 0.17), c("#440154FF", "#238A8DFF", "#FDE725FF"))

# Use the layer_fun to decide the size of the dots. 
# Within the grid.circle, we specify the radius r= sqrt(pindex(mat$percent_mat, i, j)) of the circle to be the square root of the percentage
# so the area size of the dots correspond to the percentage.

layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = "gray", fill = NA))
    grid.circle(x=x,y=y,r= sqrt(pindex(mat$percent_mat, i, j)) * unit(4, "mm"),
                gp = gpar(fill = col_fun(pindex(mat$exp_mat, i, j)), col = NA))}
  
hp<- Heatmap(mat$exp_mat,
             heatmap_legend_param=list(title= "Expression   "),
             column_title = "SRF",
             column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             width = unit(15, "cm"), # customize width of entire plot
             height = unit(5, "cm"), # customize height of entire plot
             column_names_rot = 45, # customize rotation of column names
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 12), # customize size of row label text
             border = "black",
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             row_names_side  = "left")

# Make the legend for the dot size. 
# Make sure the size is the same as the in the layer_fun which is unit(4, "mm"). You can change the size as you want.

lgd_list = list(
    Legend( labels = c(10, 25, 50, 75), title = "Percentage",
            graphics = list(
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1)  * unit(4, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(4, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(4, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(4, "mm"),
                                               gp = gpar(fill = "black"))),
            labels_gp = gpar(fontsize = 10), 
            labels_rot = 0,
            legend_height = NULL, legend_width = NULL,
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            row_gap = unit(3, "mm"),
            grid_width = unit(10, "mm")
            ))
    
# Now draw the plot
saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm"),
     annotation_legend_side = "right")
options(saved) # restore old settings

# Use more granular cell types after making some modifications

# Assuming your Seurat object is named `hecaNH`
# Create a copy of the Seurat object
hecaNH_copy <- hecaNH

# Remove cells with the 'Other' level in the copy
hecaNH_copy <- subset(hecaNH_copy, idents = "Other", invert = TRUE)

# Update the 'celltypeCollapsed' column to pool levels in the copy
hecaNH_copy$celltypeCollapsed <- as.character(hecaNH_copy$celltypeCollapsed)

# Pool 'FibF3' and 'FibF2' into one level
hecaNH_copy$celltypeCollapsed[hecaNH_copy$celltypeCollapsed %in% c('FibF3', 'FibF2')] <- 'Fibroblast'

# Pool 'SMC' and 'Periv' into one level
hecaNH_copy$celltypeCollapsed[hecaNH_copy$celltypeCollapsed %in% c('SMC', 'Periv')] <- 'SMC_Periv'

# Pool 'Endo' and 'EndoLymph' into one level
hecaNH_copy$celltypeCollapsed[hecaNH_copy$celltypeCollapsed %in% c('Endo', 'EndoLymph')] <- 'Endothelial'

# Pool 'Endo' and 'EndoLymph' into one level
hecaNH_copy$celltypeCollapsed[hecaNH_copy$celltypeCollapsed %in% c('ImmLymph', 'ImmMyel')] <- 'Immune'

# Convert back to a factor (if needed)
hecaNH_copy$celltypeCollapsed <- factor(hecaNH_copy$celltypeCollapsed)

# Update the identities based on the new cell type column in the copy
Idents(hecaNH_copy) <- hecaNH_copy$celltypeCollapsed

# Check the new levels to confirm
levels(hecaNH_copy$celltypeCollapsed)


head(hecaNH_copy)

# Let's make a cell type and phase-specific dotplot follwing instructions from https://divingintogeneticsandgenomics.com/post/how-to-make-a-multi-group-dotplot-for-single-cell-rnaseq-data/

# Customized multi-group dotplot
# We need to get the percentage of positive cells and average expression by group.
# For a single gene, put the groups into multiple rows, and each column is a cell type.

# group1 is the cell type/cluster annotation 
# group2 is any condition you want to further group, in this case, the phase

#This first part is just copy and paste from the link to prepare the function; no customization

GetMatrixFromSeuratByGroupSingle<- function(obj, feature, group1, group2){
  if (length(feature) != 1){
          stop("please only provide only one gene name")
  }
  exp_mat<- obj@assays$RNA$data[feature, ,drop=FALSE]
  count_mat<- obj@assays$RNA$counts[feature,,drop=FALSE ]
  
  meta<- obj@meta.data %>%
  tibble::rownames_to_column(var = "cell")
        
  # get the average expression matrix
  exp_df<- as.matrix(exp_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta) %>%
    group_by(gene,{{group1}}, {{group2}}) %>%
    summarise(average_expression = mean(expression)) %>%
    tidyr::pivot_wider(names_from = {{group1}}, 
                       values_from= average_expression) 
  
  exp_mat<- exp_df[, -c(1,2)] %>% as.matrix()
  rownames(exp_mat)<- exp_df %>% pull({{group2}})
  
  # get the percentage positive cell matrix
  count_df<- as.matrix(count_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    left_join(meta) %>%
    group_by(gene, {{group1}}, {{group2}}) %>%
    summarise(percentage = mean(count >0)) %>%
    tidyr::pivot_wider(names_from = {{group1}}, 
                       values_from= percentage) 

  percent_mat<- count_df[, -c(1,2)] %>% as.matrix()
  rownames(percent_mat)<- count_df %>% pull({{group2}})
  
  if (!identical(dim(exp_mat), dim(percent_mat))) {
    stop("the dimension of the two matrice should be the same!")
  }
  
  if(! all.equal(colnames(exp_mat), colnames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  
  if(! all.equal(rownames(exp_mat), rownames(percent_mat))) {
    stop("column names of the two matrice should be the same!")
  }
  return(list(exp_mat = exp_mat, percent_mat = percent_mat))
}



# Let’s get the matrices for one gene

# This is where we enter our object name, feature of interest, and condition names from the object metadata

mat<- GetMatrixFromSeuratByGroupSingle(obj= hecaNH_copy, 
                                 feature = "SRF", 
                                 celltypeCollapsed,
                                 Binary_Stage)

# take a look at the matrices

# 1. the average expression for each cell type per condition
mat$exp_mat

# 2. the percentage of cells positive for the gene of interest for each cell type per condition
mat$percent_mat

# Now, Let’s visualize it using ComplexHeatmap
# Always explore the data range before you decide how to map the data values to colors.

quantile(mat$exp_mat, c(0.1, 0.5, 0.8, 0.9))

# In this case, 0 will be mapped to #FDE725FF, 0.145 will be mapped to #238A8DFF and 0.21887 will be mapped to #440154FF. The value in-between will be linearlly interpolated to get corresponding colors

col_fun<-  circlize::colorRamp2(c(0, 0.145, 0.21887), c("#440154FF", "#238A8DFF", "#FDE725FF"))

# Use the layer_fun to decide the size of the dots. 
# Within the grid.circle, we specify the radius r= sqrt(pindex(mat$percent_mat, i, j)) of the circle to be the square root of the percentage
# so the area size of the dots correspond to the percentage.

layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = "gray", fill = NA))
    grid.circle(x=x,y=y,r= sqrt(pindex(mat$percent_mat, i, j)) * unit(4, "mm"),
                gp = gpar(fill = col_fun(pindex(mat$exp_mat, i, j)), col = NA))}
  
hp<- Heatmap(mat$exp_mat,
             heatmap_legend_param=list(title= "Expression   "),
             column_title = "SRF",
             column_title_gp = gpar(fontsize = 12, fontface = "bold"),
             width = unit(15, "cm"), # customize width of entire plot
             height = unit(5, "cm"), # customize height of entire plot
             column_names_rot = 45, # customize rotation of column names
             col=col_fun,
             rect_gp = gpar(type = "none"),
             layer_fun = layer_fun,
             row_names_gp = gpar(fontsize = 12), # customize size of row label text
             border = "black",
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             row_names_side  = "left")

# Make the legend for the dot size. 
# Make sure the size is the same as the in the layer_fun which is unit(4, "mm"). You can change the size as you want.

lgd_list = list(
    Legend( labels = c(10, 25, 50, 75), title = "Percentage",
            graphics = list(
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.1)  * unit(4, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.25) * unit(4, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.5) * unit(4, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = sqrt(0.75) * unit(4, "mm"),
                                               gp = gpar(fill = "black"))),
            labels_gp = gpar(fontsize = 10), 
            labels_rot = 0,
            legend_height = NULL, legend_width = NULL,
            title_gp = gpar(fontsize = 10, fontface = "bold"),
            row_gap = unit(3, "mm"),
            grid_width = unit(10, "mm")
            ))
    
# Now draw the plot
saved <- options(repr.plot.width = 12, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
draw(hp, annotation_legend_list = lgd_list, ht_gap = unit(1, "cm"),
     annotation_legend_side = "right")
options(saved) # restore old settings

# Make bar plot by extracting the data per instructions found here: https://bioinformatics.stackexchange.com/questions/19531/bar-graph-of-expression-data-from-seurat-object
data.df <- hecaNH@meta.data[, c("celltypeCollapsed", "Binary_Stage")] # extracts data including cell type and phase columns
data.df["value"] <- hecaNH$RNA$data["SRF", ] # adds a column for expression value of gene of interest
data.df <- data.df %>% mutate(group = paste(lineage, Binary_Stage, sep=":")) # creases a new column called "group" with both cell type and phase information

# Optional: remove cells that do not express the gene
#data.df <- data.df %>% filter(value > 0)

# create new df with mean and sd values of the above data
data.df.2 <- data.df %>% group_by(group) %>% summarise(sd=sd(value),value=mean(value)) 

# split the "group" column data back into cell type and phase columns
X <- str_split(data.df.2$'group', ":", 2, simplify=TRUE) 
data.df.2$'lineage' <- X[, 1]
data.df.2$'Binary_Stage' <- X[, 2]

#reorder factors
data.df.2$lineage <- factor(data.df.2$lineage,levels=c("Epithelial", "Mesenchymal", "Endothelial", "Immune"))
data.df.2$Binary_Stage <- factor(data.df.2$Binary_Stage,levels=c("Proliferative", "Secretory", "Menstrual"))

# Create the bar plot
saved <- options(repr.plot.width = 18, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
ggplot(data.df.2, aes(fill=Binary_Stage, y=value, x=lineage)) + 
    geom_bar(color="black", position="dodge", stat="identity")
options(saved) # restore old settings
    
#If you want to add sd error bars, add the following string after the final parenthesis above (without the quotation marks): "+ geom_errorbar(aes(ymin=value, ymax=value+sd), width=.2, position=position_dodge(.9))"

# Make bar plot by extracting the data per instructions found here: https://bioinformatics.stackexchange.com/questions/19531/bar-graph-of-expression-data-from-seurat-object
data.df <- hecaNH@meta.data[, c("lineage", "Binary_Stage")] # extracts data including cell type and phase columns
data.df["value"] <- hecaNH$RNA$data["SRF", ] # adds a column for expression value of gene of interest
data.df <- data.df %>% mutate(group = paste(lineage, Binary_Stage, sep=":")) # creases a new column called "group" with both cell type and phase information

# Optional: remove cells that do not express the gene
#data.df <- data.df %>% filter(value > 0)

# create new df with mean and sd values of the above data
data.df.2 <- data.df %>% group_by(group) %>% summarise(sd=sd(value),value=mean(value)) 

# split the "group" column data back into cell type and phase columns
X <- str_split(data.df.2$'group', ":", 2, simplify=TRUE) 
data.df.2$'lineage' <- X[, 1]
data.df.2$'Binary_Stage' <- X[, 2]

#reorder factors
data.df.2$lineage <- factor(data.df.2$lineage,levels=c("Epithelial", "Mesenchymal", "Endothelial", "Immune"))
data.df.2$Binary_Stage <- factor(data.df.2$Binary_Stage,levels=c("Proliferative", "Secretory", "Menstrual"))

# Create the bar plot
saved <- options(repr.plot.width = 18, repr.plot.height = 6) # Make the plots bigger from here on out, save old options
ggplot(data.df.2, aes(fill=Binary_Stage, y=value, x=lineage)) + 
    geom_bar(color="black", position="dodge", stat="identity") + geom_errorbar(aes(ymin=value, ymax=value+sd), width=.2, position=position_dodge(.9))
options(saved) # restore old settings
    


# Read in the processed RDS file from Sara Grimm, which contains the GD3.5 mouse (SRF flox and SRF KO) seurat object for projection to the human (HECA non-hormone treated).
# See prep code from Sara Grimm: "step4a-project_mouse_to_human-hecaNH.Rtxt"


objM2H <- readRDS("/ddn/gs1/home/marquardtrm/PRcreSRF_scRNAseq_Sara/data_rds/RM-adult_SrfKO-mouse_to_hecaNH_projection.02apr2025.rds")

objM2H
head(objM2H)


DimPlot(objM2H, reduction="ref.umap", group.by="CellType", label=T, label.size=2, repel=TRUE, raster=TRUE)+NoLegend();

levels(objM2H)

table(objM2H$genotype,objM2H$CellType);
