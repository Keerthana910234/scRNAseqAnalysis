#Integrating different RNAseq dataset which were treated with different conditions. Importing required packages or installing them.
# [Harmony](https://portals.broadinstitute.org/harmony/articles/quickstart.html) is one of the integration methods we want to use for the analysis.

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
if (!(py_module_available("scanorama"))){
  message("Installing 'Scanorama' python module...")
  py_install("scanorama", pip = T)
}
scanorama <- import('scanorama')
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
#Inputting raw 10X data for each sample before QC
#####
# Main

# Define paths
pathNaive1 <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_077_naive/filtered_feature_bc_matrix"
pathNaive2 <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_078_naive/filtered_feature_bc_matrix"
pathSotA <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_078_SotA/filtered_feature_bc_matrix" 
pathGem <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/count_SA2_078_Gem/filtered_feature_bc_matrix"

# Initialize Seurat objects
naive1 <- initSeuratObject(pathNaive1, "Naive-1")
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
createAndSaveQCPlots(naive1, qcPlotFilePath, "naive1", 10, 7000, 2000)
createAndSaveQCPlots(naive2, qcPlotFilePath, "naive2", 10, 7000, 2000)
createAndSaveQCPlots(SotA, qcPlotFilePath, "SotA", 8, 7000, 1500)
createAndSaveQCPlots(Gem, qcPlotFilePath, "Gem", 20, 7000, 2000)
#OR plot all the samples based on the parameters
plotAll(params)

naive1Filtered <- applyFilter(naive1, params$naive1)
naive2Filtered <- applyFilter(naive2, params$naive2)
SotAFiltered <- applyFilter(SotA, params$SotA)
GemFiltered <- applyFilter(Gem, params$Gem)

naive1PreProcessed <- preprocessSeurat(naive1Filtered, n_pcs = 50)
ElbowPlot(naive1PreProcessed, ndims = 50)

naive2PreProcessed <- preprocessSeurat(naive2Filtered, n_pcs = 50)
ElbowPlot(naive2PreProcessed, ndims = 50)

sotAPreProcessed <- preprocessSeurat(SotAFiltered, n_pcs = 50)
ElbowPlot(sotAPreProcessed, ndims = 50)

GemPreProcessed <- preprocessSeurat(GemFiltered, n_pcs = 50)
ElbowPlot(GemPreProcessed, ndims = 50)

combinedSeuratObject <- merge(naive2PreProcessed, y = c(sotAPreProcessed, GemPreProcessed), add.cell.ids=c("naive", "sotA", "gem"), project = "Merged",  merge.data = TRUE)
combinedSeuratObjectPreprocessed <- preprocessSeurat(combinedSeuratObject, n_pcs = 50)
ElbowPlot(combinedSeuratObjectPreprocessed, ndims = 50)

p1 <- DimPlot(object = combinedSeuratObjectPreprocessed, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = combinedSeuratObjectPreprocessed, features = "PC_1", group.by = "orig.ident", pt.size = .1)
c <- p1 + p2
c
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/beforeIntegration.png",c, height = 8, width = 12)
preProcessedSamplePath <- "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/preProcessedSingleSample"
saveRDS(naive2PreProcessed, paste0(preProcessedSamplePath, "naive2.rds"))
saveRDS(GemPreProcessed, paste0(preProcessedSamplePath, "gem.rds"))
saveRDS(sotAPreProcessed, paste0(preProcessedSamplePath, "sot.rds"))
combinedSeuratObjectPreprocessed <- FindNeighbors(combinedSeuratObjectPreprocessed, dims = 1:50, reduction = "pca")
combinedSeuratObjectPreprocessed <- FindClusters(combinedSeuratObjectPreprocessed, resolution = 0.5)
combinedSeuratObjectPreprocessed <- RunUMAP(combinedSeuratObjectPreprocessed, dims = 1:50, reduction = "pca")
saveRDS(combinedSeuratObjectPreprocessed, paste0(preProcessedSamplePath, "combined.rds"))

#analysing data before integration
#############
#Analysing data before integration
#Combining RNA-sq data with barcode data
combinedSeuratObjectPreprocessed <- readRDS( paste0(preProcessedSamplePath, "combined.rds"))
singletCells <- read_csv("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/singletList.csv")
colnames(singletCells) <- c("cellID", "barcode") 
collapsedData <- aggregate(barcode ~ cellID, data = singletCells, FUN = function(x) paste(x, collapse = ", "))

seuratCells <- data.frame(cellID = Cells(combinedSeuratObjectPreprocessed))
commonCellIDList <- merge(seuratCells, collapsedData, by = "cellID")

barcodeSeuratObjectBefore <- subset(combinedSeuratObjectPreprocessed, cells = commonCellIDList$cellID)
subsetSingletCells <- collapsedData[collapsedData$cellID %in% commonCellIDList$cellID, ]
rownames(subsetSingletCells) <- subsetSingletCells$cellID
barcodeSeuratObjectBefore$barcode <- subsetSingletCells$barcode

#Colouring clones within each sample

barcodeSampleAllBefore <- data.frame(
  sample = barcodeSeuratObjectBefore$orig.ident,
  barcode = barcodeSeuratObjectBefore$barcode,
  cellID = rownames(barcodeSeuratObjectBefore@meta.data)  # Assuming rownames are the cell IDs
) %>%
  group_by(sample, barcode) %>%
  # Calculate count but keep each cell ID
  mutate(Count = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  # Assign unique ranks within each sample based on frequency, breaking ties by their order of appearance
  mutate(clusterRank = dense_rank(-Count)) %>%
  ungroup() 

# max_cluster_index_sample1 <- max(barcodeSampleAllBefore$clusterRank[barcodeSampleAllBefore$sample == "Naive-2"])
# barcodeSampleAll <- barcodeSampleAll %>%
#   mutate(clusterRank = ifelse(sample == "Sot", clusterRank + max_cluster_index_sample1, clusterRank))
# max_cluster_index_sample2 <- max(barcodeSampleAll$clusterRank[barcodeSampleAll$sample == "Sot"])
# barcodeSampleAll <- barcodeSampleAll %>%
#   mutate(clusterRank = ifelse(sample == "Gem", clusterRank + max_cluster_index_sample2, clusterRank))
barcodeSampleAllBefore <-  barcodeSampleAllBefore %>%
  column_to_rownames(var = "cellID")

barcodeRankedSeuratObjectBefore <- barcodeSeuratObjectBefore
barcodeRankedSeuratObjectBefore$cloneRank <- barcodeSampleAllBefore$clusterRank
umapBefore <- FetchData(barcodeSeuratObjectBefore, vars = c("umap_1", "umap_2", "barcode", "orig.ident"))
umapBefore$cellID <- rownames(umapBefore)
plotLargestBarcodeBeforeIntegration <- ggplot(umapBefore, aes(x = umap_1, y = umap_2)) +
  geom_point(data = ~ subset(., barcode != "AGCATGCCCGCCAGCGCTACTCCGCACCTGCCACCCTCAGCTCTCGGGCC"), color = "lightgrey", size = 2, alpha = 0.4) +
  geom_point(data = ~ subset(., barcode == "AGCATGCCCGCCAGCGCTACTCCGCACCTGCCACCCTCAGCTCTCGGGCC"), color = "navyblue", size = 2, alpha = 1) +
  # geom_text_repel(data = subset(umapDataAllBarcodes, allCondition == 1), aes(label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  labs(title = paste("UMAP of cells present in all 3 conditions")) +
  theme_void() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), # Center the title
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
   facet_wrap(~orig.ident)
plotLargestBarcodeBeforeIntegration

##########################################################################################################################################
##########################################################################################################################################
library(Seurat)
library(ggplot2)
samples <- unique(barcodeSeuratObjectBefore$orig.ident)

# Retrieve unique sample identifiers
samples <- unique(barcodeSeuratObjectBefore$orig.ident)

# Define a vibrant color palette for the highlighted barcodes; adjust the number based on unique barcodes
barcodeColors <- grDevices::rainbow(length(unique(barcodeRankedSeuratObjectBefore$barcode)))

# Name the colors by barcode for consistent mapping
names(barcodeColors) <- unique(barcodeRankedSeuratObjectBefore$barcode)

# Light grey for background cells
barcodeColors <- c("lightgrey", barcodeColors)

# Iterate over each sample
# Iterate over each sample
# Iterate over each sample
# Define a color palette avoiding light colors
paletteName <- "Dark2"  # Change this as needed

for (sample in samples) {
  for (size in 1:3) {
    barcodeRankedSeuratObjectBefore$highlight <- ifelse(
      barcodeRankedSeuratObjectBefore$orig.ident == sample & barcodeRankedSeuratObjectBefore$cloneRank == size,
      barcodeRankedSeuratObjectBefore$barcode,
      "Other"
    )
    
    uniqueBarcodes <- unique(barcodeRankedSeuratObjectBefore$barcode[barcodeRankedSeuratObjectBefore$highlight != "Other"])
    numberColors <- min(length(uniqueBarcodes), brewer.pal.info[paletteName, "maxcolors"])
    barcodeColors <- setNames(brewer.pal(numberColors, paletteName), uniqueBarcodes)
    barcodeColors <- c("lightgrey", barcodeColors)  # Add light grey for 'Other'
    
    umap_data <- FetchData(barcodeRankedSeuratObjectBefore, vars = c("umap_1", "umap_2", "highlight"))
    
    plot <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = highlight)) +
      geom_point(data = ~ subset(., highlight == "Other"), color = "lightgrey", size = 2, alpha = 0.4) +
      geom_point(data = ~ subset(., highlight != "Other"), aes(color = highlight), size = 2, alpha = 1) +
      scale_color_manual(values = barcodeColors) +
      labs(title = paste("UMAP for Sample", sample, "CloneRank", size)) +
      theme_minimal() +
      theme(legend.position = "none") +
      theme_void() +
      theme(legend.position = "none", 
            plot.background = element_blank(), # This sets the plot background to be blank
            panel.background = element_blank(), # Ensure the panel background is blank
            panel.border = element_blank(), # Removes panel border if present
            panel.grid.major = element_blank(), # Removes major grid lines
            panel.grid.minor = element_blank()) # Removes minor grid lines # Adjust alpha for point transparency if desired
    
    
    print(plot)
    ggsave(plot = plot, filename = paste0("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/beforeIntegration/cloneWise_", sample, "_", size, ".png"))
  }
}
plotAll <- DimPlot(barcodeRankedSeuratObjectBefore, group.by = "orig.ident")
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/beforeIntegration/samplePlot.png", plotAll)
plotAll
########################################################################################################################################################################################
# Integrating with Harmony
########################################################################################################################################################################################
options(repr.plot.height = 2.5, repr.plot.width = 6)
listFeatures <- Features(integratedObject)
integratedObject <- RunHarmony(combinedSeuratObjectPreprocessed, c("orig.ident"), plot_convergence = T)
integratedObject <- JoinLayers(integratedObject)

p1after <- DimPlot(object = integratedObject, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2after <- VlnPlot(object = integratedObject, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
d <- p1after + p2after
d
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/afterIntegration.png",d, height = 8, width = 12)

integratedObject <- FindNeighbors(integratedObject, dims = 1:50, reduction = "harmony")
integratedObject <- FindClusters(integratedObject, resolution = 0.5)
integratedObject <- RunUMAP(integratedObject, dims = 1:50, reduction = "harmony")
left <- DimPlot(integratedObject, reduction = "umap")
right <- DimPlot(integratedObject, reduction = "umap", group.by = "orig.ident")
e <- left + right
e
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/clusterUMAP.png",e, height = 8, width = 12)
saveRDS(integratedObject, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/integratedDataHarmony.rds")
markersDE <- FindAllMarkers(integratedObject)
#Analysing integrated data
##########################################################################################################################################
integratedObject <- readRDS("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/integratedDataHarmony.rds")

#Combining RNA-sq data with barcode data
singletCells <- read_csv("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/singletList.csv")
colnames(singletCells) <- c("cellID", "barcode") 
collapsedData <- aggregate(barcode ~ cellID, data = singletCells, FUN = function(x) paste(x, collapse = ", "))

seuratCells <- data.frame(cellID = Cells(integratedObject))
commonCellIDList <- merge(seuratCells, collapsedData, by = "cellID")

barcodeSeuratObject <- subset(integratedObject, cells = commonCellIDList$cellID)
subsetSingletCells <- collapsedData[collapsedData$cellID %in% commonCellIDList$cellID, ]
rownames(subsetSingletCells) <- subsetSingletCells$cellID
barcodeSeuratObject$barcode <- subsetSingletCells$barcode

saveRDS(barcodeSeuratObject, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/barcodeSeuratObject.rds")
###################################################################################################################
# Plotting common clones between conditions
###################################################################################################################
barcodeSeuratObject <- readRDS("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/barcodeSeuratObject.rds")
barcodeSamples <- data.frame(
  barcode = barcodeSeuratObject$barcode,
  sample = barcodeSeuratObject$orig.ident
)

# Define sample pairs
samplePairs <- list(
  "Naive-2 and Gem" = c("Naive-2", "Gem"),
  "Gem and Sot" = c("Gem", "Sot"),
  "Naive-2 and Sot" = c("Naive-2", "Sot")
)

pairwiseIntersections <- list()  # Initialize the list to store results

# Calculate pairwise intersections
for (pairName in names(samplePairs)) {
  pair <- samplePairs[[pairName]]
  pairwiseIntersections[[pairName]] <- barcodeSamples %>%
    group_by(barcode) %>%
    filter(sample %in% pair) %>%  # Ensure barcode appears in the specified pair samples
    summarise(count = n_distinct(sample)) %>%
    filter(count == length(pair)) %>%  # Ensure barcode appears in exactly those two samples
    pull(barcode)
}


#Calculating the barcodes present in all three samples
barcodeSeuratObjectAllCondition <- barcodeSeuratObject
barcodeSamplesAllCondition <- barcodeSamples %>% 
  group_by(barcode) %>% 
  summarise(count = n_distinct(sample)) %>%
 filter(count ==3)

#Plotting only the barcodes that are present in all 3
barcodeSeuratObjectAllCondition$allCondition = 0
barcodeSeuratObjectAllCondition$allCondition[barcodeSeuratObjectAllCondition$barcode %in% barcodeSamplesAllCondition$barcode] <- 1
umapDataAllBarcodes <- FetchData(barcodeSeuratObjectAllCondition, vars = c("umap_1", "umap_2", "allCondition", "barcode", "orig.ident"))
umapDataAllBarcodes$cellID <- rownames(umapDataAllBarcodes)
umapData_all_conditions <- umapDataAllBarcodes[umapDataAllBarcodes$allCondition == 1, ]

# Assign numbers to the barcodes present in all conditions
umapData_all_conditions$BarcodeNumber <- as.numeric(factor(umapData_all_conditions$barcode))

# Now combine the data with the filtered BarcodeNumbers with the original data
umapDataAllBarcodes <- merge(umapDataAllBarcodes, umapData_all_conditions[, c("BarcodeNumber", "cellID")], by = "cellID", all.x = T)
# umapDataAllBarcodes$BarcodeNumber <- as.numeric(factor(umapDataAllBarcodes$barcode))

plotColouredBySample <- ggplot(umapDataAllBarcodes, aes(x = umap_1, y = umap_2)) +
  geom_point(data = ~ subset(., allCondition == 0), color = "lightgrey", size = 2, alpha = 0.4) +
  geom_point(data = ~ subset(., allCondition == 1), aes(color = orig.ident), size = 2, alpha = 1) +
  geom_text_repel(data = subset(umapDataAllBarcodes, allCondition == 1), aes(label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  labs(title = paste("UMAP of cells present in all 3 conditions")) +
  theme_void() +
  theme(legend.position = "right", 
        plot.title = element_text(hjust = 0.5), # Center the title
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + facet_wrap(~ orig.ident) 

plotColouredBySample

#####PLotting in PCA Space
pcaDataAllBarcodes <- FetchData(barcodeSeuratObjectAllCondition, vars = c("harmony_1", "harmony_2", "allCondition", "barcode", "orig.ident"))
pcaDataAllBarcodes$cellID <- rownames(pcaDataAllBarcodes)
pcaData_all_conditions <- pcaDataAllBarcodes[pcaDataAllBarcodes$allCondition == 1, ]

# Assign numbers to the barcodes present in all conditions
pcaData_all_conditions$BarcodeNumber <- as.numeric(factor(pcaData_all_conditions$barcode))

# Now combine the data with the filtered BarcodeNumbers with the original data
pcaDataAllBarcodes <- merge(pcaDataAllBarcodes, pcaData_all_conditions[, c("BarcodeNumber", "cellID")], by = "cellID", all.x = T)
# umapDataAllBarcodes$BarcodeNumber <- as.numeric(factor(umapDataAllBarcodes$barcode))

plotColouredBySamplePCA <- ggplot(pcaDataAllBarcodes, aes(x = harmony_1, y = harmony_2)) +
  geom_point(data = ~ subset(., allCondition == 0), color = "lightgrey", size = 2, alpha = 0.4) +
  geom_point(data = ~ subset(., allCondition == 1), aes(color = orig.ident), size = 2, alpha = 1) +
  # geom_text_repel(data = subset(pcaDataAllBarcodes, allCondition == 1), aes(label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  labs(title = paste("PCA 1 and PCA2 cells present in all 3 conditions")) +
  # theme_void() +
  # theme(legend.position = "right", 
  #       plot.title = element_text(hjust = 0.5), # Center the title
  #       plot.background = element_blank(), 
  #       panel.background = element_blank(), 
  #       panel.border = element_blank(), 
  #       panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank()) +
  facet_wrap(~ orig.ident) 

plotColouredBySamplePCA
ggsave(plotColouredBySamplePCA, filename ="/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/harmonyUMAP/pcaPlotsAllConditions.png")

plotColouredByBarcode <- ggplot(umapDataAllBarcodes, aes(x = umap_1, y = umap_2)) +
  geom_point(data = ~ subset(., allCondition == 0), color = "lightgrey", size = 2, alpha = 0.4) +
  geom_point(data = ~ subset(., allCondition == 1), aes(color = barcode), size = 2, alpha = 1) +
  geom_text_repel(data = subset(umapDataAllBarcodes, allCondition == 1), aes(label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  labs(title = paste("UMAP of cells present in all 3 conditions")) +
  theme_void() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), # Center the title
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plotColouredByBarcode


plotColouredByBarcodeNoLabel <- ggplot(umapDataAllBarcodes, aes(x = umap_1, y = umap_2)) +
  geom_point(data = ~ subset(., allCondition == 0), color = "lightgrey", size = 2, alpha = 0.4) +
  geom_point(data = ~ subset(., allCondition == 1), aes(color = barcode), size = 2, alpha = 1) +
  # geom_text_repel(data = subset(umapDataAllBarcodes, allCondition == 1), aes(label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  labs(title = paste("UMAP of cells present in all 3 conditions")) +
  theme_void() +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5), # Center the title
        plot.background = element_blank(), 
        panel.background = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

plotColouredByBarcodeNoLabel

background_data <- umapDataAllBarcodes %>%
  # Create a copy of each point for each unique orig.ident to ensure it appears in all facets
  crossing(origIdent = unique(umapDataAllBarcodes$orig.ident)) %>%
  mutate(orig.ident = origIdent) %>%
  select(-origIdent)

colors <- c("#8dd3c7",
            "#bc80bd",
            "#ffffb3",
            "blue",
            "#fb8072",
            "darkgreen",
            "gold",
            "#b3de69",
            "hotpink",
            "red")
plotColouredByBarcodeNoLabelFacet <- ggplot() +
  # Add grey background points for all data
  geom_point(data = background_data, aes(x = umap_1, y = umap_2), color = "lightgrey", size = 2, alpha = 0.2) +
  # Add colored points only where condition is met
  geom_point(data = umapDataAllBarcodes %>% filter(allCondition == 1),
             aes(x = umap_1, y = umap_2, color = barcode), size = 2, alpha = 1) +
  # geom_text_repel(data = subset(umapDataAllBarcodes, allCondition == 1), aes(x = umap_1, y = umap_2, label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  scale_color_manual(values = colors) +
  facet_wrap(~ orig.ident) +
  labs(title = "UMAP of cells present in all 3 conditions") +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +  
  theme(aspect.ratio = 1/1.2)

plotColouredByBarcodeNoLabelFacet

ggsave(plot = plotColouredByBarcodeNoLabelFacet, filename = "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/harmonyUMAP/commonClonesAcrossConditions.png", width = 15, height = 6)
combinedPlot <- cowplot::plot_grid(plotColouredBySample, plotColouredByBarcode, plotColouredByBarcodeNoLabel, nrow = 1)
combinedPlot
write_csv(umapDataAllBarcodes, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plotData/umapDataAllBarcodesCommonBarcodes.csv")

umapDataAllBarcodesCount <- umapDataAllBarcodes %>%
  group_by(barcode, orig.ident) %>%  # Group by barcode and sample
  mutate(sizePerBarcodePerSample = n()) %>%  # Count number of entries per group
  ungroup()  #

facet_labels <- umapDataAllBarcodesCount %>%
  group_by(orig.ident) %>%
  summarise(title = paste("Counts per Barcode: ", paste(BarcodeNumber, count, sep=": ", collapse=", "), sep=""))

temp = umapDataAllBarcodesCount %>% filter(allCondition ==1)
unique_barcodes <- unique(temp$BarcodeNumber)

# Map colors to these barcodes
# Note: This assumes the number of unique barcodes does not exceed the length of your colors list
color_mapping <- setNames(colors, unique_barcodes)

# Add the color for the background data
color_mapping["background"] = "lightgrey"

plotColouredByBarcodeSizeNoLabelFacet <- function(number, color) {
  # Filter the data for the specific BarcodeNumber
  data_to_plot <- umapDataAllBarcodesCount %>%
    filter(BarcodeNumber == number) %>%
    mutate(BarcodeNumber = as.factor(BarcodeNumber),  # Convert BarcodeNumber to factor
           facet_label = paste(orig.ident, "\nSize of clone: ", sizePerBarcodePerSample))
  
  plot <- ggplot(data_to_plot, aes(x = umap_1, y = umap_2)) +
    geom_point(data = background_data, color = "lightgrey", size = 2, alpha = 0.2) +
    # geom_point(data = filter(data_to_plot, allCondition == 1),
    #            aes(color = BarcodeNumber), size = 2, alpha = 1) +
    geom_point(data = filter(data_to_plot, allCondition == 1),
              aes(color = orig.ident), size = 2, alpha = 1) +
    # scale_color_manual(values = color_mapping) +
    # facet_wrap(~ facet_label, scales = "free") +
    labs(title = paste("UMAP of Clone Number ", number, " Present in All 3 Conditions")) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5),
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1/1.2
    )
  savePath = "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/umapCommonClones/"
  ggsave(plot = plot, filename = paste0(savePath, "sampleWiseCloneNumber", number, ".png") )
  print(plot)
}

colorsNew <- c("green3",
            "purple",
            "navyblue",
            "blue",
            "#fb8072",
            "darkgreen",
            "gold",
            "violet",
            "hotpink3",
            "red")

# Call the function with the corrected color mapping
plotColouredByBarcodeSizeNoLabelFacet(1, color_mapping)
for(number in 1:10){
plotColouredByBarcodeSizeNoLabelFacet(number, colorsNew[number]) 
  
}
plot

for (pairName in names(samplePairs)) {
  pair <- samplePairs[[pairName]]
  pairwiseIntersections[[pairName]] <- barcodeSamples %>%
    group_by(barcode) %>%
    summarise(count = n_distinct(sample))
}
# Debug print to check the results
print(pairwiseIntersections)
barcodeSamples <- barcodeSamples %>%
  rownames_to_column(var = "RowName")
barcodeSamples <- barcodeSamples %>%
  mutate(clusterIndex = as.integer(factor(barcode, levels = unique(barcode))))

barcodeSeuratObject$barcodeID <- barcodeSamples$clusterIndex
library(ggrepel)
plotSubset <- function(seuratObject, barcodes, title, pair) {
  # Determine valid cells based on barcodes and pair criteria
  validCells <- which(seuratObject@meta.data$barcode %in% barcodes & seuratObject@meta.data$orig.ident %in% pair)
  
  # Assign a highlight column in the metadata to differentiate valid from non-valid cells
  seuratObject$highlight <- ifelse(seq_along(seuratObject@meta.data$barcode) %in% validCells, seuratObject@meta.data$barcode[validCells], "Other")
  
  # Extract UMAP data for all cells
  umap_data <- FetchData(seuratObject, vars = c("umap_1", "umap_2", "barcode", "highlight", "barcodeID"))
  
  # Define colors using hcl.colors for better granularity and more options
  uniqueBarcodes <- unique(umap_data$barcode[umap_data$highlight != "Other"])
  barcodeColors <- hcl.colors(length(uniqueBarcodes), "Dark 3")  # +1 for the 'Other' category
  barcodeColors <- c("lightgrey", barcodeColors)
  names(barcodeColors) <- c("Other", uniqueBarcodes)
  
  # Create the plot
  plot <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = highlight)) +
    geom_point(data = umap_data[umap_data$highlight == "Other", ], aes(color = "Other"), size = 2, alpha = 0.4) +
    geom_point(data = umap_data[umap_data$highlight != "Other", ], aes(color = barcode), size = 2, alpha = 1) +
    geom_text_repel(data = umap_data[umap_data$highlight != "Other", ], aes(label = barcodeID, color = barcode), size = 3, box.padding = 0.35, point.padding = 0.5) +
    scale_color_manual(values = barcodeColors) +
    ggtitle(title) +
    theme_minimal() +
    theme(legend.position = "none")
  
  print(plot)
}

# Execute plotting for each pair
for (pairName in names(pairwiseIntersections)) {
  barcodes <- pairwiseIntersections[[pairName]]
  plotSubset(barcodeSeuratObject, barcodes, sprintf("Intersection of %s", pairName), samplePairs[[pairName]])
}


##########
#Colouring clones within each sample
barcodeSampleAll <- data.frame(
  sample = barcodeSeuratObject$orig.ident,
  barcode = barcodeSeuratObject$barcode,
  cellID = rownames(barcodeSeuratObject@meta.data)  # Assuming rownames are the cell IDs
) %>%
  group_by(sample, barcode) %>%
  # Calculate count but keep each cell ID
  mutate(Count = n()) %>%
  ungroup() %>%
  group_by(sample) %>%
  # Assign unique ranks within each sample based on frequency, breaking ties by their order of appearance
  mutate(clusterRank = dense_rank(-Count)) %>%
  ungroup() 
# Printing the final resu

max_cluster_index_sample1 <- max(barcodeSampleAll$clusterRank[barcodeSampleAll$sample == "Naive-2"])
barcodeSampleAll <- barcodeSampleAll %>%
  mutate(clusterRank = ifelse(sample == "Sot", clusterRank + max_cluster_index_sample1, clusterRank))
max_cluster_index_sample2 <- max(barcodeSampleAll$clusterRank[barcodeSampleAll$sample == "Sot"])
barcodeSampleAll <- barcodeSampleAll %>%
  mutate(clusterRank = ifelse(sample == "Gem", clusterRank + max_cluster_index_sample2, clusterRank))
barcodeSampleAll <-  barcodeSampleAll %>%
  column_to_rownames(var = "cellID")
barcodeRankedSeuratObject <- barcodeSeuratObject
barcodeRankedSeuratObject$cloneRank <- barcodeSampleAll$clusterRank

samples <- unique(barcodeRankedSeuratObject$orig.ident)

# Retrieve unique sample identifiers
samples <- unique(barcodeRankedSeuratObject$orig.ident)

# Define a vibrant color palette for the highlighted barcodes; adjust the number based on unique barcodes
barcodeColors <- grDevices::rainbow(length(unique(barcodeRankedSeuratObject$barcode)))

# Name the colors by barcode for consistent mapping
names(barcodeColors) <- unique(barcodeRankedSeuratObject$barcode)

# Light grey for background cells
barcodeColors <- c("lightgrey", barcodeColors)

# Iterate over each sample
# Iterate over each sample
# Iterate over each sample
# Define a color palette avoiding light colors
paletteName <- "Dark2"  # Change this as needed

for (sample in samples) {
  for (size in 1:3) {
    barcodeRankedSeuratObject$highlight <- ifelse(
      barcodeRankedSeuratObject$orig.ident == sample & barcodeRankedSeuratObject$cloneRank == size,
      barcodeRankedSeuratObject$barcode,
      "Other"
    )
    
    uniqueBarcodes <- unique(barcodeRankedSeuratObject$barcode[barcodeRankedSeuratObject$highlight != "Other"])
    numberColors <- min(length(uniqueBarcodes), brewer.pal.info[paletteName, "maxcolors"])
    barcodeColors <- setNames(brewer.pal(numberColors, paletteName), uniqueBarcodes)
    barcodeColors <- c("lightgrey", barcodeColors)  # Add light grey for 'Other'
    
    umap_data <- FetchData(barcodeRankedSeuratObject, vars = c("umap_1", "umap_2", "highlight"))
    
    plot <- ggplot(umap_data, aes(x = umap_1, y = umap_2, color = highlight)) +
      geom_point(data = ~ subset(., highlight == "Other"), color = "lightgrey", size = 2, alpha = 0.4) +
      geom_point(data = ~ subset(., highlight != "Other"), aes(color = highlight), size = 2, alpha = 1) +
      scale_color_manual(values = barcodeColors) +
      labs(title = paste("UMAP for Sample", sample, "CloneRank", size)) +
      theme_minimal() +
      theme(legend.position = "none") +
      theme_void() +
      theme(legend.position = "none", 
            plot.background = element_blank(), # This sets the plot background to be blank
            panel.background = element_blank(), # Ensure the panel background is blank
            panel.border = element_blank(), # Removes panel border if present
            panel.grid.major = element_blank(), # Removes major grid lines
            panel.grid.minor = element_blank()) # Removes minor grid lines # Adjust alpha for point transparency if desired
    
    
    print(plot)
    ggsave(plot = plot, filename = paste0("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/harmonyUMAP/cloneWise_", sample, "_", size, ".png"))
  }
}


######
#Plotting individual clones in different conditions



facet_labels <- umapDataAllBarcodesAllCount %>%
  group_by(orig.ident) %>%
  summarise(title = paste("CloneSize: ", paste(BarcodeNumber, sizePerBarcodePerSample, sep=": ", collapse=", "), sep=""))

background_data <- umapDataAllBarcodesAllCount %>%
  # Create a copy of each point for each unique orig.ident to ensure it appears in all facets
  crossing(barcodeNumber = unique(umapDataAllBarcodesAllCount$BarcodeNumber)) %>%
  mutate(BarcodeNumber = barcodeNumber) %>%
  select(-barcodeNumber)

umapDataAllBarcodesAllCount <- umapDataAllBarcodesCount
umapDataAllBarcodesAllCount$BarcodeNumber <- as.numeric(factor(umapDataAllBarcodesAllCount$barcode))
umapDataAllBarcodesAllCount$facet_label <- facet_labels


umapDataAllBarcodesAllCount <- umapDataAllBarcodesAllCount %>%
  mutate(BarcodeNumber = as.factor(BarcodeNumber),  # Convert BarcodeNumber to factor
         facet_label = paste(BarcodeNumber, "\nSize of clone: ", sizePerBarcodePerSample))

plotSampleWiseCloneWise <- ggplot() +
  # Add grey background points for all data
  # Add colored points only where condition is met
  geom_point(data = umapDataAllBarcodesAllCount %>% filter(orig.ident == "Sot")%>% filter(sizePerBarcodePerSample >2), aes(x = umap_1, y = umap_2, color = barcode), size = 2, alpha = 1) +
  geom_point(data = umapDataAllBarcodesAllCount %>%
               filter(orig.ident != "Sot" | sizePerBarcodePerSample <= 2),
             aes(x = umap_1, y = umap_2), color = "lightgrey", size = 2, alpha = 0.2) +
  # geom_text_repel(data = subset(umapDataAllBarcodes, allCondition == 1), aes(x = umap_1, y = umap_2, label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(title = "UMAP of clones in Sot") +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +  
  theme(aspect.ratio = 1/1.2)
# plotSampleWiseCloneWise <- plotSampleWiseCloneWise +   geom_point(data = umapDataAllBarcodesAllCount, aes(x = umap_1, y = umap_2), color = "lightgrey", size = 2, alpha = 0.2) 
plotSampleWiseCloneWise


plotSampleWiseCloneWise <- ggplot() +
  # Add grey background points for all data
  # Add colored points only where condition is met
  geom_point(data = umapDataAllBarcodesAllCount %>% filter(orig.ident == "Sot")%>% filter(sizePerBarcodePerSample >2), aes(x = umap_1, y = umap_2, color = barcode), size = 2, alpha = 1) +
  geom_point(data = umapDataAllBarcodesAllCount %>%
               filter(orig.ident != "Sot" | sizePerBarcodePerSample <= 2),
             aes(x = umap_1, y = umap_2), color = "lightgrey", size = 2, alpha = 0.2) +
  # geom_text_repel(data = subset(umapDataAllBarcodes, allCondition == 1), aes(x = umap_1, y = umap_2, label = BarcodeNumber), size = 4, box.padding = 0.35, point.padding = 0.5) +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(title = "UMAP of clones in Sot") +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +  
  theme(aspect.ratio = 1/1.2)
# plotSampleWiseCloneWise <- plotSampleWiseCloneWise +   geom_point(data = umapDataAllBarcodesAllCount, aes(x = umap_1, y = umap_2), color = "lightgrey", size = 2, alpha = 0.2) 
plotSampleWiseCloneWise

ggsave(plotSampleWiseCloneWise, filename = "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/harmonyUMAP/clonesSotBackground.png", height = 5, width = 5)
#################################################################################################################################################################
#################################################################################################################################################################
#Scanorama integration
group_names <- c("Naive-2", "Gem", "Sot")

datasets <- list()
gene_lists <- list()

for (group_name in group_names) {
  # Get the expression data for the specified assay
  assay_data <- combinedSeuratObjectPreprocessed[["RNA"]][group_name]
  print(dim(assay_data))
  # Transpose the data matrix (convert to matrix if not already)
  transposed_data <- t(as.matrix(assay_data))
  print(dim(transposed_data))
  group_indices <- which(Idents(combinedSeuratObjectPreprocessed) == group_name)
  # Store results in lists
  datasets[[group_name]] <- transposed_data
  gene_lists[[group_name]] <- rownames(assay_data)
}

names(datasets) = NULL
names(gene_lists) = NULL

integrated_corrected_data = scanorama$correct(datasets, gene_lists, return_dimred = TRUE, return_dense = TRUE, ds_names = group_names, verbose = TRUE)

corrected_scanorama <- t(do.call(rbind, integrated_corrected_data[[2]]))
colnames(corrected_scanorama) <- colnames(combinedSeuratObjectPreprocessed)
rownames(corrected_scanorama) <- integrated_corrected_data[[3]]
corrected_scanorama_pca <- t(do.call(rbind, integrated_corrected_data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(combinedSeuratObjectPreprocessed)

scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
combinedSeuratObjectPreprocessed[["scanorama"]] <- scanorama_assay
DefaultAssay(combinedSeuratObjectPreprocessed) <- "scanorama"



all.genes <- rownames(combinedSeuratObjectPreprocessed)
scanoramaSeuratObject <- ScaleData(combinedSeuratObjectPreprocessed, features = all.genes)
scanoramaSeuratObject <- RunPCA(object = scanoramaSeuratObject, assay = "scanorama", features= all.genes, reduction.name = "pca_scanorama")

scanoramaSeuratObject <- FindNeighbors(object=scanoramaSeuratObject, features = all.genes, dims = 1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
scanoramaSeuratObject <- FindClusters(object=scanoramaSeuratObject,graph.name = "scanorama_snn", resolution = 0.5)
scanoramaSeuratObject <- RunUMAP(object = scanoramaSeuratObject, reduction = "pca_scanorama", dims = 1:50, reduction.name = "umap_scanorama") 
plotScanorama <- DimPlot(scanoramaSeuratObject, reduction = "pca_scanorama",group.by = "orig.ident", pt.size = 1) +ggtitle("after scanorama")
plotScanorama
plotBefore <- DimPlot(combinedSeuratObjectPreprocessed, group.by = "orig.ident", pt.size = 1) +ggtitle("before integration")
plot <- plotScanorama + plotBefore
plot
saveRDS(scanoramaSeuratObject, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/integratedDataScanorama.rds")
scanoramaSeuratObject <- readRDS("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/integratedData/integratedDataScanorama.rds")
pcaScanorama <- DimPlot(scanoramaSeuratObject,reduction = "pca_scanorama", group.by = "orig.ident")
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/scanoramaPCA.png", pcaScanorama)
umap <- DimPlot(scanoramaSeuratObject, reduction = "umap", group.by = "orig.ident")
umap
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/umap.png", umapScanorama)
