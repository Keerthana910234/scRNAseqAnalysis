
#Created by Keerthana M Arun on 20240618 and last modified by Keerthana M Arun on 20240618 at 3:27 PM
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = F)
# Load necessary libraries
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis) 
library(RColorBrewer)
library(reticulate)
library(readr)
library(ggrepel)
library(SeuratWrappers)
library(scCustomize)
py_install("scvi")
py_install("scanpy")
scvi <- import("scvi")
# Function definitions
##############################################################
# Function to filter Seurat objects using QC parameters
filterQC <- function(seuratObject, minCounts = 0, minGenes = 0, maxMitoPercentage = 100, maxCounts = Inf, maxGenes = Inf) {
  seuratObject <- subset(seuratObject, subset = nFeature_RNA > minGenes & nFeature_RNA < maxGenes &
                           percent.mito < maxMitoPercentage & nCount_RNA > minCounts & nCount_RNA < maxCounts)
  return(seuratObject)
}

#Function to plot QC plots for all samples in one go
plotAll <- function(params) {
  lapply(params, function(p) {
    createAndSaveQCPlots(p$seuratObject, qcPlotFilePath, p$prefix, 
                         p$percentMitoCutoff, p$minCounts, p$minFeatures,
                         p$maxCounts, p$maxFeatures)
  })
}

# Function to create and save QC plots
createAndSaveQCPlots <- function(seuratObject, path, prefix, percentageMitoCutoff, minimumCounts, minimumFeatures, maximumCounts = NA, maximumFeatures = NA) {
  # Scatter plot of RNA counts vs. mitochondrial percentage
  plot1 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "percent.mito") + 
    geom_vline(xintercept = minimumCounts, linetype = "dotted", color = "red") + 
    geom_hline(yintercept = percentageMitoCutoff, linetype = "dashed", color = "blue")
  if (!is.na(maximumCounts)) {
    plot1 <- plot1 + geom_vline(xintercept = maximumCounts, linetype = "dashed", color = "red")
  }
  
  # Scatter plot of RNA features vs. mitochondrial percentage
  plot2 <- FeatureScatter(seuratObject, feature1 = "nFeature_RNA", feature2 = "percent.mito") + 
    geom_vline(xintercept = minimumFeatures, linetype = "dotted", color = "red") + 
    geom_hline(yintercept = percentageMitoCutoff, linetype = "dashed", color = "blue")
  if (!is.na(maximumFeatures)) {
    plot2 <- plot2 + geom_vline(xintercept = maximumFeatures, linetype = "dashed", color = "red")
  }
  
  # Scatter plot of RNA counts vs. RNA features
  plot3 <- FeatureScatter(seuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
    geom_vline(xintercept = minimumCounts, linetype = "dotted", color = "red") + 
    geom_hline(yintercept = minimumFeatures, linetype = "dotted", color = "red")
  if (!is.na(maximumCounts) && !is.na(maximumFeatures)) {
    plot3 <- plot3 + geom_vline(xintercept = maximumCounts, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = maximumFeatures, linetype = "dashed", color = "red")
  }
  
  # Violin plot for mitochondrial gene percentage
  plot4 <- VlnPlot(seuratObject, features = c("percent.mito")) + 
    geom_hline(yintercept = percentageMitoCutoff, linetype = "dashed", color = "blue")
  
  # Save plots to specified path
  ggsave(paste0(path, prefix, "CountVsPercentMito.png"), plot = plot1)
  ggsave(paste0(path, prefix, "FeatureVsPercentMito.png"), plot = plot2)
  ggsave(paste0(path, prefix, "CountVsFeature.png"), plot = plot3)
  ggsave(paste0(path, prefix, "PercentMitoViolin.png"), plot = plot4)
}

# Function to read and initialize Seurat objects from 10X data
initSeuratObject <- function(path, projectName) {
  counts <- Read10X(path)
  seuratObject <- CreateSeuratObject(counts, min.cells = 3, min.features = 200, project = projectName)
  seuratObject[["percent.mito"]] <- PercentageFeatureSet(seuratObject, pattern = "^MT-")
  return(seuratObject)
}

# Function to apply filtering for a specific dataset
applyFilter <- function(seuratObject, params) {
  filterQC(seuratObject, minCounts = params$minCounts, minGenes = params$minFeatures, 
           maxMitoPercentage = params$percentMitoCutoff, maxCounts = params$maxCounts, maxGenes = params$maxFeatures)
}

#Function to perform pre-processing on the cells
preprocessSeurat <- function(seuratObject, n_pcs = 30){
  seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
  seuratObject <- ScaleData(seuratObject, features = Features(seuratObject))
  seuratObject <- RunPCA(seuratObject, n_pcs = n_pcs, features = Features(seuratObject))
  return(seuratObject)
}

plotUmap <- function(dataframe, xAxis, yAxis, colourBy, legend = "right"){
  # Convert string inputs into symbols and then into expressions
  p <- ggplot(data = dataframe, aes(x = !!rlang::sym(xAxis), 
                                    y = !!rlang::sym(yAxis), 
                                    color = !!rlang::sym(colourBy))) +
    geom_point(alpha = 0.8, size = 1) +  # Adding point layer with some transparency
    theme_void() +
    theme(legend.position = legend,  # Adjust this to put legend where you prefer
          plot.title = element_text(hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 0.6)
  return(p)
}

plotPca <- function(dataframe, xAxis, yAxis, colourBy, legend = "right"){
  # Convert string inputs into symbols and then into expressions
  p <- ggplot(data = dataframe, aes(x = !!rlang::sym(xAxis), 
                                    y = !!rlang::sym(yAxis), 
                                    color = !!rlang::sym(colourBy))) +
    geom_point(alpha = 0.8, size = 1) +  # Adding point layer with some transparency
    labs( x = xAxis,  # Label for x-axis
         y = yAxis) +  # Label for y-axis
    theme_minimal() +  # Using minimal theme which includes axes
    theme(legend.position = legend,  # Adjust this to put legend where you prefer
          plot.title = element_text(hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_rect(fill = "white", colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          panel.grid.major = element_line(linewidth = 0.1, linetype = 'solid', colour = "gray"),
          panel.grid.minor = element_blank(),
          aspect.ratio = 0.6)
  return(p)
}

library(reticulate)
devtools::install_github("cellgeni/sceasy")
library(loomR)
library(sceasy)
library(Seurat)
library(patchwork)
h5ad_file <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/scviIntegrated_20240621.h5ad"
sceasy::convertFormat(
  h5ad_file,
  from = "anndata",
  to = "seurat",
  outFile = "~/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/scviFullGene.rds",
  main_layer = "counts"
)
scvi <- readRDS('~/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/scviFullGene.rds')
scviMaybe[["RNA"]] <- as(object = scviMaybe[["RNA"]], Class = "Assay5")

scviMaybe <- FindNeighbors(scviMaybe, reduction = "scVI", dims = 1:30)
scviMaybe <- FindClusters(scviMaybe, resolution = 0.5, cluster.name = "scvi_clusters")
scviMaybe <- RunUMAP(scviMaybe, dims = 1:30, reduction = "scVI")

umapData <- FetchData(scviMaybe, vars = c("umap_1", "umap_2", "batch", "scvi_clusters"))
umapData$cellID <- rownames(umapData)

plotUmapSample <- plotUmap(umapData, "umap_1", "umap_2", "batch", legend = "bottom")
plotUmapSample
ggsave("~/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsUMAP/scVISamples.png", plotUmapSample)
write_csv(umapData,"~/Keerthana/subiaHanxiaoDataAnalysis/plotData/integrationPlotsUMAP/scvi_0.5resolution.csv")
plotUmapSample <- plotUmap(umapData, "umap_1", "umap_2", "scvi_clusters", legend = "none")
plotUmapSample
ggsave("~/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsUMAP/scVICluster_0.5resolution.png", plotUmapSample)
saveRDS(scviMaybe, "~/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/scVI_R.rds")

pcaData <- FetchData(scviMaybe, vars = c("SCVI_1", "SCVI_2", "batch", "scvi_clusters"))
plotPCASample <- plotPca(pcaData, "SCVI_1", "SCVI_2", "batch", legend = "bottom")
plotPCASample
ggsave("~/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsPCA/scVISample.png", plotPCASample)
pcaData$cellID <- rownames(pcaData)
write_csv(pcaData,"~/Keerthana/subiaHanxiaoDataAnalysis/plotData/integrationPlotsPCA/scvi.csv")

#Build the rds object for scVI integrated h5ad object from parts
library(data.table)
allPath <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/rdsObjects/scVIcomponents/"
scVIDirectory <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/rdsObjects/scVIcomponents/"
expr_matrix <- fread(paste0(allPath, "matrix.csv"), data.table = F, header = TRUE)
expr_matrix <- expr_matrix[, -1]


metadata <- read.csv(paste0(allPath, 'metadata.csv'))
features <- read.csv(paste0(allPath, 'features.csv'), row.names = 1)
scvi <- read.csv(paste0(allPath, 'scviEmbedding.csv'), row.names = 1)
rownames(scvi) <- rownames(metadata) 
# Use the data to create a Seurat object if needed
library(Seurat)
seurat_object <- CreateSeuratObject(counts = t(as.matrix(expr_matrix)), project = "scvi")
rownames(seurat_object) <- rownames(features)
colnames(seurat_object) <- rownames(metadata)
seurat_object <- AddMetaData(seurat_object, metadata = metadata)

colnames(scvi) <- paste("scVI", 1:ncol(scvi), sep="_")

# Adding to Seurat object
seurat_object[['scvi']] <- CreateDimReducObject(
  embeddings = as.matrix(scvi),
  key = "scVI_",
  assay = "RNA"
)

seurat_object <- FindNeighbors(seurat_object, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution = 0.5, cluster.name = "scvi_clusters")
seurat_object <- RunUMAP(seurat_object, reduction = "scvi", dims = 1:10)

umapData <- FetchData(seurat_object, vars = c("umap_1", "umap_2", "batch", "scvi_clusters"))
umapData$cellID <- rownames(umapData)
write_csv(umapData, )
plotUmapSample <- plotUmap(umapData, "umap_1", "umap_2", "batch", legend = "bottom")
plotUmapSample
write_csv(umapData,"~/Keerthana/subiaHanxiaoDataAnalysis/plotData/integrationPlotsUMAP/scvi_new.csv")

ggsave("~/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsUMAP/scVISample_0.5resolution_new.png", plotUmapSample)
pcaData <- FetchData(seurat_object, vars = c("scVI_1", "scVI_2", "batch", "scvi_clusters"))
pcaData$cellID <- rownames(pcaData)
plotPCASample <- plotPca(pcaData, "scVI_1", "scVI_2", "batch", legend = "bottom")
plotPCASample
ggsave("~/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsPCA/scVISample_new.png", plotPCASample)
saveRDS(seurat_object, "~/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/scVI_new.rds")

write_csv(pcaData,"~/Keerthana/subiaHanxiaoDataAnalysis/plotData/integrationPlotsPCA/scvi_new.csv")
