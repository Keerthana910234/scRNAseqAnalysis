# [Harmony](https://portals.broadinstitute.org/harmony/articles/quickstart.html) is one of the integration methods we want to use for the analysis.

#Created by Keerthana M Arun on 20240617 and last modified by Keerthana M Arun on 20240619 at 4:00 PM

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
    geom_point(alpha = 0.8) +  # Adding point layer with some transparency
    labs(title ="UMAP of samples integrated by Harmony") +
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
    geom_point(alpha = 0.8) +  # Adding point layer with some transparency
    labs(title ="PCA of samples integrated by Harmony",  # Adjusted title for PCA
         x = xAxis,  # Label for x-axis
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

#####
# Main

# Define paths
# pathNaive1 <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_077_naive/filtered_feature_bc_matrix"
pathNaive2 <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_078_naive/filtered_feature_bc_matrix"
pathSotA <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_078_SotA/filtered_feature_bc_matrix" 
pathGem <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_078_Gem/filtered_feature_bc_matrix"

# Initialize Seurat objects
# naive1 <- initSeuratObject(pathNaive1, "Naive-1")
naive2 <- initSeuratObject(pathNaive2, "Naive-2")
SotA <- initSeuratObject(pathSotA, "Sot")
Gem <- initSeuratObject(pathGem, "Gem")

# Define QC plot file path
qcPlotFilePath <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/qcPlots/"

params <- list(
  # naive1 = list(seuratObject = naive1, prefix = "naive1", percentMitoCutoff = 10, minCounts = 7000, minFeatures = 2000, maxCounts = Inf, maxFeatures = Inf),
  naive2 = list(seuratObject = naive2, prefix = "naive2", percentMitoCutoff = 10, minCounts = 7000, minFeatures = 2000, maxCounts = Inf, maxFeatures = Inf),
  SotA = list(seuratObject = SotA, prefix = "SotA", percentMitoCutoff = 8, minCounts = 7000, minFeatures = 1500, maxCounts = Inf, maxFeatures = Inf),
  Gem = list(seuratObject = Gem, prefix = "Gem", percentMitoCutoff = 20, minCounts = 7000, minFeatures = 2000, maxCounts = Inf, maxFeatures = Inf)
)

# Create and save plots for both samples
# createAndSaveQCPlots(naive1, qcPlotFilePath, "naive1", 10, 7000, 2000)
# createAndSaveQCPlots(naive2, qcPlotFilePath, "naive2", 10, 7000, 2000)
# createAndSaveQCPlots(SotA, qcPlotFilePath, "SotA", 8, 7000, 1500)
# createAndSaveQCPlots(Gem, qcPlotFilePath, "Gem", 20, 7000, 2000)
#OR plot all the samples based on the parameters
plotAll(params)

naive2Filtered <- applyFilter(naive2, params$naive2)
SotAFiltered <- applyFilter(SotA, params$SotA)
GemFiltered <- applyFilter(Gem, params$Gem)

combinedSeuratObject <- merge(naive2Filtered, y = c(SotAFiltered, GemFiltered), add.cell.ids=c("naive", "sot", "gem"), project = "Merged",  merge.data = TRUE)
combinedSeuratObjectPreprocessed <- preprocessSeurat(combinedSeuratObject, n_pcs = 50)
ElbowPlot(combinedSeuratObjectPreprocessed, ndims = 50)

f = c("MT-ND1")
# , "AC078923.1", , "MT-ATP6", "MT-ND2", "MT-CYB", "BCAR3", "MT-ATP8", "MT-ND4", "CD44", 
# "HPCAL1")
FeaturePlot(combinedSeuratObjectPreprocessed, features = f)

# integratedObject <- RunHarmony(combinedSeuratObjectPreprocessed, c("orig.ident"), plot_convergence = T)
integratedObject <- IntegrateLayers(object = combinedSeuratObjectPreprocessed, method = HarmonyIntegration, orig.reduction = "pca",
                                          new.reduction = 'harmony', features = Features(combinedSeuratObjectPreprocessed))

integratedObject <- JoinLayers(integratedObject)

listFeatures <- Features(integratedObject)

p1after <- DimPlot(object = integratedObject, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
# p2after <- VlnPlot(object = integratedObject, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
p1after
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/afterIntegration.png",d, height = 8, width = 12)

integratedObject <- FindNeighbors(integratedObject, dims = 1:50, reduction = "harmony")
integratedObject <- FindClusters(integratedObject, resolution = 0.5)
integratedObject <- RunUMAP(integratedObject, dims = 1:50, reduction = "harmony")

umapData <- FetchData(integratedObject, vars = c("umap_1", "umap_2", "orig.ident", "seurat_clusters"))
umapData$cellID <- rownames(umapData)

pcaData <- FetchData(integratedObject, vars = c("harmony_1", "harmony_2", "orig.ident", "seurat_clusters"))

plotUmapCluster <- plotUmap(umapData, "umap_1", "umap_2", "seurat_clusters", legend = "none")
plotUmapCluster
write_csv(umapData, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plotData/integrationPlotsUMAP/harmonyClustersNew_0.5resolution.csv")
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsUMAP/harmonyClustersNew_0.5resolution.png",plotUmapCluster, height = 8, width = 12)
plotUmapSample <- plotUmap(umapData, "umap_1", "umap_2", "orig.ident", legend = "bottom") 
plotUmapSample
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsUMAP/harmonySamplesNew.png",plotUmapSample, height = 8, width = 12)

plotPcaSample <- plotPca(pcaData, "harmony_1", "harmony_2", "orig.ident")
plotPcaSample
plotPcaCluster <- plotPca(pcaData, "harmony_1", "harmony_2", "seurat_clusters")
plotPcaCluster
write_csv(pcaData, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plotData/integrationPlotsPCA/harmonyClustersNew_0.5resolution.csv")

ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsPCA/harmonySamplesNew.png",plotPcaSample, height = 8, width = 12)
saveRDS(integratedObject,"/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/rdsObjects/harmonyNew.rds")
integratedObject <- readRDS("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/rdsObjects/harmonyNew.rds")
