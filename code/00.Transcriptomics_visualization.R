
### Sevda Rafatov
## 01.08.2023
## Transcriptomics project visualization

###
#We conducted a transcriptomics project using the Seurat package in R to analyze single-cell RNA sequencing data. 
#The data were visualized to explore the cellular composition and gene expression patterns in different cell types. 
#We used various plotting techniques and palettes to enhance the clarity and aesthetics of the visualizations.
#The scRNA reference data was annotated with cell types, and we renamed the clusters to reflect the corresponding cell identities. 
#To facilitate the interpretation of the results, we utilized the "umap" reduction method for dimensionality 
#reduction and applied the "polychrome" palette with 21 distinct colors to represent different cell types.
#Initially, we obtained a UMAP plot showing all cells with labeled cell types. Next, we created a grouped UMAP plot based on the original identities of the cells. 
#The legend levels were reordered to represent three different experimental conditions, namely "hu140_8.6 pcw," "hu084_9.0 pcw," and "hu122_10.7 pcw."
#To further investigate gene expression patterns in specific cell types, we curated a list of genes of interest. 
#We verified the presence of these genes in the data and proceeded to create violin plots for each gene using the "polychrome" palette. 
#These plots provide valuable insights into the distribution and expression levels of the selected genes across different cell types.
#Moreover, we generated single-cell heatmaps for the chosen gene set. The heatmap visualizations showcase the expression patterns of the selected genes in individual cells, 
#providing a detailed view of gene activity within different cell types.
###

# Load required libraries
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(forcats)

# Importing pre-processed data in the Seurat object
my_data <- readRDS("/mnt/DATA_4TB/projects/transcriptomics_Heather/new_results/snRNA_tissue_integ4(maxGfilter)AfterClust.rds")

# Renaming the clusters with the corresponding cell types
new.cluster.ids <- c("CMV", "EnTh", "CmVCnD", "CmA4", "CmG2M", "Fb1", "EndoCd", "Fb3","CmA1", "FbVIC",  "CmA3",    "Fb2", "SMC", "CmA2a", "EpiCd", "CmA2b", "CmMit", "Bl",  "SCP", "NE", "EnVL")
names(new.cluster.ids) <- levels(my_data)
my_data <- RenameIdents(my_data, new.cluster.ids)

# Display information about the annotated cell types
cat("scRNA reference data annotated with clusters cell types\n")

# Create a new column 'celltype' in the Seurat object to store the cell type annotations
my_data[["celltype"]] <- my_data@active.ident

# Set the default assay to "integrated" to work with the integrated dataset
DefaultAssay(my_data) <- "integrated"

# Obtain a UMAP plot of the cells with polychrome palette for color representation
DimPlot(my_data, reduction = "umap", label = TRUE, repel = FALSE, cols = DiscretePalette(n=21, palette = "polychrome"))

# Obtain a grouped UMAP plot based on the original identities of the cells
p1 <- DimPlot(my_data, reduction = "umap", group.by = "orig.ident", cols = DiscretePalette(n=21, palette = "polychrome"))

# Reorder levels in the legend to represent different experimental conditions
p1$data$orig.ident <- fct_relevel(p1$data$orig.ident, "hu140_8.6 pcw", "hu084_9.0 pcw", "hu122_10.7 pcw") 
p1

# Create a combined plot with the grouped UMAP and the individual UMAP side by side
p2 <- DimPlot(my_data, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

# Visualize each condition side by side using the split.by argument in DimPlot
p3 <- DimPlot(my_data, reduction = "umap", split.by = "orig.ident", cols = DiscretePalette(n=21, palette = "polychrome"))

# Reorder levels in the legend for the split.by plot
p3$data$orig.ident <- fct_relevel(p3$data$orig.ident, "hu140_8.6 pcw", "hu084_9.0 pcw", "hu122_10.7 pcw") 

# Display the split.by plot with reordered legend levels
p3

#Changin original ident to hu140_8.6 pcw

# Assuming you have already loaded and processed your data into 'my_data'

# Assuming you have already loaded and processed your data into 'my_data'

# Specify the orig.ident level you want to select
selected_orig_ident <- "hu140_8.6 pcw"

# Find the cell indices that match the desired orig.ident level
selected_indices <- which(my_data$orig.ident == selected_orig_ident)

# Create a new Seurat object containing only the selected cells
selected_seurat <- my_data[selected_indices, ]

# Now, 'selected_seurat' contains only the cells from the desired orig.ident level

#Change the assay to RNA
DefaultAssay(my_data) <- "RNA"  
# Create a gene vector with genes of interest
genes <- c('PRRX1', 'EDIL3', 'FABP4', 'MYL7','NRF2','PSTN','SDF1','ALDH1A1', 'ALDH1A2', 'ASTN2', 'BMP10', 'BMP4', 'CDH19', 'COL2A1', 'COL6A6', 'CRABP1', 'CRABP2', 'CXCL12', 'CXCL14', 'CXCR4', 'DLK1', 'EBF2', 'EDIL3', 'EDN1', 'ELN', 'ENPEP', 'ERBB3', 'FABP4', 'FBLN1', 'FLI1', 'FRZB', 'FST', 'GATA4', 'GATA5', 'GPC3', 'HCN4', 'HES1', 'HGF', 'HHIP', 'HOXA5', 'IGF1', 'IRX1', 'IRX2', 'IRX3', 'IRX4', 'IRX5', 'IRX6', 'ISL1', 'MDFI', 'MYH6', 'MYH7', 'MYL2', 'MYL7', 'MYOM2', 'NKX2-5', 'NNAT', 'NR2F2', 'NRG1', 'NRP2', 'NRXN1', 'PECAM1', 'PHOX2B', 'PITX2', 'PLP1', 'POSTN', 'PRRX1', 'SCN7A', 'SHOX2', 'SNAI1', 'SNAI2', 'SOX10', 'SOX4', 'SPP1', 'SPRY1', 'SST', 'SSTR2A', 'TBX1', 'TBX18', 'TBX3', 'TBX5', 'TCF21', 'TF', 'TRH', 'TWIST', 'VSNL1', 'WNT2', 'XKR4', 'ZEB2')

# Check if all genes are present in the data
missing_genes <- genes[!genes %in% rownames(my_data)]
if (length(missing_genes) > 0) {
  cat("Warning: The following genes are not found in the data:\n")
  cat(paste(missing_genes, collapse = ", "), "\n")
}

pdf(file = "multipleViolinPlot.pdf")
for (gene in genes) {
  if (gene %in% rownames(my_data)) {
    print(VlnPlot(my_data, features = gene, cols = DiscretePalette(n=21, palette = "polychrome")))
  } else {
    cat("Skipping gene", gene, "as it is not found in the data.\n")
  }
}
dev.off()   # important line to indicate the process to stop write to the file

# Create gene vector
genes <- c('PRRX1', 'EDIL3', 'FABP4', 'MYL7','NRF2','PSTN','SDF1','ALDH1A1', 'ALDH1A2', 'ASTN2', 'BMP10', 'BMP4', 'CDH19', 'COL2A1', 'COL6A6', 'CRABP1', 'CRABP2', 'CXCL12', 'CXCL14', 'CXCR4', 'DLK1', 'EBF2', 'EDIL3', 'EDN1', 'ELN', 'ENPEP', 'ERBB3', 'FABP4', 'FBLN1', 'FLI1', 'FRZB', 'FST', 'GATA4', 'GATA5', 'GPC3', 'HCN4', 'HES1', 'HGF', 'HHIP', 'HOXA5', 'IGF1', 'IRX1', 'IRX2', 'IRX3', 'IRX4', 'IRX5', 'IRX6', 'ISL1', 'MDFI', 'MYH6', 'MYH7', 'MYL2', 'MYL7', 'MYOM2', 'NKX2-5', 'NNAT', 'NR2F2', 'NRG1', 'NRP2', 'NRXN1', 'PECAM1', 'PHOX2B', 'PITX2', 'PLP1', 'POSTN', 'PRRX1', 'SCN7A', 'SHOX2', 'SNAI1', 'SNAI2', 'SOX10', 'SOX4', 'SPP1', 'SPRY1', 'SST', 'SSTR2A', 'TBX1', 'TBX18', 'TBX3', 'TBX5', 'TCF21', 'TF', 'TRH', 'TWIST', 'VSNL1', 'WNT2', 'XKR4', 'ZEB2')

# Verify if all are in your data (you should have everything TRUE)
missing_genes <- genes[!genes %in% rownames(my_data)]
if (length(missing_genes) > 0) {
  cat("Warning: The following genes are not found in the data:\n")
  cat(paste(missing_genes, collapse = ", "), "\n")
}

pdf(file = "multipleFeaturePlot.pdf")
for (gene in genes) {
  if (gene %in% rownames(my_data)) {
    print(FeaturePlot(my_data, features = gene))
  } else {
    cat("Skipping gene", gene, "as it is not found in the data.\n")
  }
}
dev.off()   # important line to indicate the process to stop write to the file


# Assuming you have loaded your Seurat object named 'my_data' and manually re-clustered your cells.
# The 'new.cluster.ids' vector should contain the new cluster names you provided.

# Rename the clusters with the provided cell type labels
new.cluster.ids <- c("I", "II", "I", "I", "I", "III", "II", "III","I", "III",  "I", "III", "III", "I", "IV", "I", "I", "V",  "VI", "VI", "II")
names(new.cluster.ids) <- levels(my_data)
my_data <- RenameIdents(my_data, new.cluster.ids)

# Curate a list of genes of interest for Heatmap visualization
features <- c('BMP10', 'MYH6', 'MYOM2', 'TBX5', 'SHOX2', 'EDN1', 'POSTN', 'PECAM1', 'NRG1', 'FLI1', 'EBF2', 'EDIL3', 'COL6A6', 'ELN', 'ENPEP', 'GPC3', 'HGF', 'IGF1', 'PRRX1', 'SCN7A', 'ALDH1A1', 'ASTN2', 'CDH19', 'ISL1', 'NRXN1', 'XKR4', 'ALDH1A2', 'HHIP', 'TBX18', 'NRP2', 'MYH7', 'SPP1')

# Create a UMAP plot with 'polychrome' palette to visualize the cells
DimPlot(my_data, reduction = "umap", label = TRUE, repel = FALSE, cols = DiscretePalette(n=21, palette = "polychrome"))

# Create a DotPlot for the selected genes to visualize their expression across cell types
# The size of the dot corresponds to the percentage of cells expressing the gene in each cluster,
# and the color represents the average expression level
dot_plot <- DotPlot(my_data, features = features)

# Rotate the x-axis labels for better visualization
dot_plot + RotatedAxis() + theme(axis.text.x = element_text(size = 7))

# Create a single-cell Heatmap of gene expression for the selected features
# Down sample the data for faster computation
DoHeatmap(subset(my_data, downsample = 100), features = features, size = 3)

