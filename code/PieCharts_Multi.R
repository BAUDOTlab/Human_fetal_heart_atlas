# Contributed by Camille Humbert

#!/usr/bin/env Rscript
library(Seurat)
library(SPOTlight)
library(ggplot2)
library(scCustomize)

################
# INIT & PATHS #
################
workingDir <- "/media/gbim/Volume12TB/Volume12TB/2020-43-2021-20/RCTD_deconvolution"
RDir <- file.path(workingDir, "R_Multi")
RDataDir <- file.path(RDir,"Rdata")

spatialPath <- file.path(RDataDir, "Coeur_integration_14_06_2024.rds")
load(spatialPath)
spatial <- Coeur.integration

samplesSpatial <- c("2021A90","2021A91", "2021A17","2021A18","2021A19","2021A20")
spatialSplit<- readRDS(file.path(RDataDir, paste0("SpatialSplit_Multi_", paste(samplesSpatial, collapse = "_"), ".rds")))
ref <- readRDS(file.path(RDataDir, "snRNA_3DataIntegrated and labeled(all).rds"))
removeCluster <- "16"


cell_types <- levels(ref)
cell_types <- cell_types[!cell_types == removeCluster]
pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")

# PieCharts
system(command = paste0("mkdir -p ", file.path(RDir,"reports","PieCharts")))
for (sample in samplesSpatial){
  pos <- GetTissueCoordinates(spatialSplit[[sample]])
  colnames(pos) <- c("x","y")
  deconData <- spatialSplit[[sample]]@meta.data[,cell_types]
  pdf(file = file.path(RDir,"reports","PieCharts",paste0(sample,".pdf")),
      width = 12,
      height = 8)
  
  print(plotSpatialScatterpie(x = pos, y = deconData, pie_scale = 0.35) + scale_fill_manual(name = "",  values = pal))
  dev.off()
}
sample <- "2021A18"
pos <- GetTissueCoordinates(spatialSplit[[sample]])
colnames(pos) <- c("x","y")
deconData <- spatialSplit[[sample]]@meta.data[,cell_types]
pdf(file = file.path(RDir,"reports","PieCharts",paste0(sample,".pdf")),
    width = 12,
    height = 8)
print(plotSpatialScatterpie(x = pos, y = deconData, pie_scale = 0.55) + scale_fill_manual(name = "",  values = pal))
dev.off()


# SpatialFeaturePlot pour chaque échantillon avec la bonne résolution
sample <- "2021A91"
for (cl in cell_types){
  tempPath <- file.path(RDir,"reports","SpatialFeaturePlotDeconvolution", sample)
  system(command = paste0("mkdir -p ", tempPath))
  pdf(file = file.path(tempPath,paste0(sample,"_cluster",cl,".pdf")),
      width = 12,
      height = 8)
  print(SpatialFeaturePlot(spatialSplit[[sample]], features = cl, pt.size.factor = 2.35))
  dev.off()
}
sample <- "2021A90"
for (cl in cell_types){
  tempPath <- file.path(RDir,"reports","SpatialFeaturePlotDeconvolution", sample)
  system(command = paste0("mkdir -p ", tempPath))
  pdf(file = file.path(tempPath,paste0(sample,"_cluster",cl,".pdf")),
      width = 12,
      height = 8)
  print(SpatialFeaturePlot(spatialSplit[[sample]], features = cl, pt.size.factor = 2.35))
  dev.off()
}
sample <- "2021A17"
for (cl in cell_types){
  tempPath <- file.path(RDir,"reports","SpatialFeaturePlotDeconvolution", sample)
  system(command = paste0("mkdir -p ", tempPath))
  pdf(file = file.path(tempPath,paste0(sample,"_cluster",cl,".pdf")),
      width = 12,
      height = 8)
  print(SpatialFeaturePlot(spatialSplit[[sample]], features = cl, pt.size.factor = 2.2))
  dev.off()
}
sample <- "2021A18"
for (cl in cell_types){
  tempPath <- file.path(RDir,"reports","SpatialFeaturePlotDeconvolution", sample)
  system(command = paste0("mkdir -p ", tempPath))
  pdf(file = file.path(tempPath,paste0(sample,"_cluster",cl,".pdf")),
      width = 12,
      height = 8)
  print(SpatialFeaturePlot(spatialSplit[[sample]], features = cl, pt.size.factor = 2.6))
  dev.off()
}
sample <- "2021A19"
for (cl in cell_types){
  tempPath <- file.path(RDir,"reports","SpatialFeaturePlotDeconvolution", sample)
  system(command = paste0("mkdir -p ", tempPath))
  pdf(file = file.path(tempPath,paste0(sample,"_cluster",cl,".pdf")),
      width = 12,
      height = 8)
  print(SpatialFeaturePlot(spatialSplit[[sample]], features = cl, pt.size.factor = 1.9))
  dev.off()
}
sample <- "2021A20"
for (cl in cell_types){
  tempPath <- file.path(RDir,"reports","SpatialFeaturePlotDeconvolution", sample)
  system(command = paste0("mkdir -p ", tempPath))
  pdf(file = file.path(tempPath,paste0(sample,"_cluster",cl,".pdf")),
      width = 12,
      height = 8)
  print(SpatialFeaturePlot(spatialSplit[[sample]], features = cl, pt.size.factor = 2.05))
  dev.off()
}

# Table with average percentage of each nuclei cluster in each spatial cluster and outputs a table
system(command = paste0("mkdir -p ", file.path(RDir,"reports","Tables")))
histList <- list()
for (sample in samplesSpatial){
  clustersSpatial <- levels(spatial)
  histTable <- data.frame(matrix(nrow = length(clustersSpatial), ncol = length(cell_types)))
  colnames(histTable) <- cell_types
  rownames(histTable) <- clustersSpatial
  clustersSpatialSample <- levels(spatialSplit[[sample]])
  for (cl in clustersSpatialSample){
    sampleTemp <- spatialSplit[[sample]]
    sampleTemp <- subset(sampleTemp, subset = seurat_clusters == cl)
    histTable[cl,] <- colMeans(sampleTemp@meta.data[,cell_types], na.rm = TRUE)
  }
  histList[[sample]] <- histTable
  rownames(histTable) <- paste0("Spatial_",rownames(histTable))
  colnames(histTable) <- paste0("Nuclei_",colnames(histTable))
  histTable <- round(histTable,4)
  histTable <- tibble::rownames_to_column(histTable, "Cluster")
  write.table(histTable, file.path(RDir,"reports","Tables",paste0(sample,"_ratio.tsv")), row.names = FALSE, sep = "\t")
  histTable <- round(histTable[2:ncol(histTable)]*100, 2)
  rownames(histTable) <- paste0("Spatial_",rownames(histTable))
  histTable <- tibble::rownames_to_column(histTable, "Cluster")
  write.table(histTable, file.path(RDir,"reports","Tables",paste0(sample,"_percent.tsv")), row.names = FALSE, sep = "\t")
}


# Pie charts for averages by cluster
for (sample in samplesSpatial){
  pdf(file.path(RDir,"reports","PieCharts",paste0(sample,"_clusterAverage.pdf")),
      width = 12,
      height = 8)
  par(mfrow = c(4,4))
  par(mar = c(0,0,0,0))
  for (cl in levels(spatialSplit[[sample]])){
    temp <- as.data.frame(t(histList[[sample]][cl,])*100)
    colnames(temp) <- c("percent")
    temp <- tibble::rownames_to_column(temp, "cluster")
    temp$colours <- pal[1:nrow(temp)]
    pie(temp$percent, col = temp$colours, border = FALSE, labels = NA, main = paste0("\nCl ",cl), radius = 0.8)
  }
  dev.off()
  pdf(file.path(RDir,"reports","PieCharts",paste0(sample,"_sup5percent_clusterAverage.pdf")),
      width = 12,
      height = 8)
  par(mfrow = c(4,4))
  par(mar = c(0,0,0,0))
  for (cl in levels(spatialSplit[[sample]])){
    temp <- as.data.frame(t(histList[[sample]][cl,])*100)
    temp[temp <5] <- 0
    colnames(temp) <- c("percent")
    temp <- tibble::rownames_to_column(temp, "cluster")
    temp$colours <- pal[1:nrow(temp)]
    pie(temp$percent, col = temp$colours, border = FALSE, labels = NA, main = paste0("\nCl ",cl), radius = 0.8)
  }
  dev.off()
}
