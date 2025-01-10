# Contributed by Camille Humbert

#!/usr/bin/env Rscript
library(spacexr)
library(Matrix)
library(Seurat)
library(popkin)
library(Polychrome)
library(dplyr)
library(scCustomize)

################
# INIT & PATHS #
################
workingDir <- "/media/gbim/Volume12TB/Volume12TB/2020-43-2021-20/RCTD_deconvolution"
RDir <- file.path(workingDir, "R")
RDataDir <- file.path(RDir,"Rdata")
spatialPath <- file.path(RDataDir, "Coeur_integration_14_06_2024.rds")
mypal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
removeCluster <- "16"
samplesSpatial <- c("2021A90","2021A91","2021A17","2021A18","2021A19","2021A20")

# Integrated vs 2021A90, 2021A91, 2021A17, 2021A18, 2021A19, 2021A20
#############################################
nucleiPath <- file.path(RDataDir, "snRNA_3DataIntegrated and labeled(all).rds")
ref <- readRDS(nucleiPath)
load(spatialPath)
spatial <- Coeur.integration

#######################################
# DATA PREPROCESSING AND RUNNING RCTD #
#######################################
# Remove unwanted cluster
subsetRef <- subset(ref, subset = seurat_clusters != removeCluster)
subsetRef@meta.data$seurat_clusters <- droplevels(subsetRef@meta.data$seurat_clusters)

# Preprocessing
subsetRefNorm <- SCTransform(subsetRef)
subsetRefNorm <- RunPCA(subsetRefNorm)
subsetRefNorm <- RunUMAP(object = subsetRefNorm, dims = 1:30)
ref <- subsetRefNorm

countsRef <- ref[["RNA"]]$counts
clusterRef <- as.factor(ref$seurat_clusters)
names(clusterRef) <- colnames(ref)
nUMIRef <- ref$nCount_RNA
names(nUMIRef) <- colnames(ref)
reference <- Reference(countsRef, clusterRef, nUMIRef)
cell_types <- levels(ref)
spatialSplit <- SplitObject(spatial, split.by = "orig.ident")

for (sample in samplesSpatial){
  spatialTemp <- spatialSplit[[sample]]
  countsSp <- spatialTemp[["Spatial"]]$counts
  coordsSp <- GetTissueCoordinates(spatialTemp)
  colnames(coordsSp) <- c("x", "y")
  coordsSp[is.na(colnames(coordsSp))] <- NULL
  query <- SpatialRNA(coordsSp, countsSp, colSums(countsSp))
  RCTD <- create.RCTD(query, reference, max_cores = 8, CELL_MIN_INSTANCE = 20, MAX_MULTI_TYPES = 10)
  RCTD <- run.RCTD(RCTD, doublet_mode = "multi")
  barcodes <- rownames(RCTD@spatialRNA@coords)
  subWeightsTable <- data.frame(matrix(0L, nrow = length(barcodes), ncol = length(cell_types)))
  rownames(subWeightsTable) <- barcodes
  colnames(subWeightsTable) <- cell_types
  names(RCTD@results) <- barcodes
  for (barcode in barcodes){
    if (length(RCTD@results[[barcode]]$cell_type_list) > 1){
      for (cl in RCTD@results[[barcode]]$cell_type_list){
        subWeightsTable[barcode, cl] <- RCTD@results[[barcode]]$sub_weights[[cl]]
      }
    }else{
      cl <- RCTD@results[[barcode]]$cell_type_list
      subWeightsTable[barcode, cl] <- RCTD@results[[barcode]]$sub_weights
    }
  }
  spatialSplit[[sample]] <- AddMetaData(spatialTemp, metadata = subWeightsTable[,cell_types])
}
saveRDS(spatialSplit, file.path(RDataDir, paste0("SpatialSplit_Multi_", paste(samplesSpatial, collapse = "_"), ".rds")))

