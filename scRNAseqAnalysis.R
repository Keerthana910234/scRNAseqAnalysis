#Created by Keerthana M Arun on 20240621 and last modified by Keerthana M Arun on 20240621 at 4:52 PM


library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis) 
library(dplyr)
library(reticulate)
library(readr)
library(ggrepel)

#Functions for later use

#Plotting neat UMAPs
plotUmap <- function(dataframe, xAxis, yAxis, colourBy, legend = "right"){
  # Convert string inputs into symbols and then into expressions
  p <- ggplot(data = dataframe, aes(x = !!rlang::sym(xAxis), 
                                    y = !!rlang::sym(yAxis), 
                                    color = !!rlang::sym(colourBy))) +
    geom_point(alpha = 0.8, size = 1) +  # Adding point layer with some transparency
    # labs(title ="UMAP of samples integrated by Harmony") +
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

#Plotting neat PCA plots
plotPca <- function(dataframe, xAxis, yAxis, colourBy, legend = "right"){
  # Convert string inputs into symbols and then into expressions
  p <- ggplot(data = dataframe, aes(x = !!rlang::sym(xAxis), 
                                    y = !!rlang::sym(yAxis), 
                                    color = !!rlang::sym(colourBy))) +
    geom_point(alpha = 0.8, size = 1) +  # Adding point layer with some transparency
    # labs(title ="PCA of samples integrated by Harmony",  # Adjusted title for PCA
    #      x = xAxis,  # Label for x-axis
    #      y = yAxis) +  # Label for y-axis
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

randomTesting <- function(allGenes, deGenes, bulkDeGenes, n = 10000){
  # Set seed for reproducibility
  set.seed(42)
  # Define the number of bootstrap samples
  numSamples <- n
  # Define the size of each sample based on the number of genes in your comparison list
  sampleSize <- length(deGenes)
  # Perform bootstrap sampling
  overlapCounts <- replicate(numSamples, {
    # Sample genes without replacement
    sampledGenes <- sample(allGenes, sampleSize, replace = FALSE)
    # Calculate intersection with your gene list
    length(intersect(sampledGenes, bulkDeGenes))
  })
  # Observed number of intersected genes (make sure this is calculated or defined correctly)
  observedIntersection <- length(intersect(deGenes, bulkDeGenes))
  # Calculate p-value
  p_value <- sum(overlapCounts >= observedIntersection) / numSamples
  # Output results
  cat("Observed Intersection:", observedIntersection, "\n")
  cat("P-Value for observing equal or more extreme overlap by chance:", p_value, "\n")
  return(p_value)
}

#Input will be the integrated rds object - from Harmony
seuratObject <- readRDS("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/rdsObjects/harmonyNew.rds")

umapData <- FetchData(seuratObject, vars = c("umap_1", "umap_2", "orig.ident", "seurat_clusters"))
umapData$cellID <- rownames(umapData)


plotUmapSample <- plotUmap(umapData, "umap_1", "umap_2", "orig.ident", legend = "none") 
plotUmapSample
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/integrationPlotsUMAP/harmonySamples_plain.svg",plotUmapSample, height = 8, width = 8)



#Run DEseq on different conditions

#First subset the genes that are also present in bulk RNAseq data
bulkRNAseqGeneList <- read_csv("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/deSeqBulkRnaSeq_withoutOneGem.csv")

#Find common genes in both bulk and the rds object
commonGene <- data.frame(geneName = intersect(bulkRNAseqGeneList$gene, Features(seuratObject)))
#Subset the Features in the seuratObject to the ones in bulkRNAseqGeneList
seuratObject <- subset(seuratObject, features = commonGene$geneName)

#Differential expression testing between different conditions

Idents(seuratObject) <- "orig.ident"
sotarasibVsNaive <- FindMarkers(seuratObject, ident.1 = "Naive-2", ident.2 = "Sot", assay="RNA", fc.slot = "counts")
sotarasibVsNaive <- read_csv("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/sotarasibVsNaive_scRNAseq.csv")
sotarasibVsNaive <- sotarasibVsNaive %>% filter(p_val_adj < 0.001, abs(avg_log2FC) >= 2.5)
sotarasibVsNaive$gene <- rownames(sotarasibVsNaive)
write_csv(sotarasibVsNaive, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/sotarasibVsNaive_scRNAseq.csv")

gemcitabineVsNaive <- FindMarkers(seuratObject, ident.1 = "Naive-2", ident.2 = "Gem")
gemcitabineVsNaive <- read_csv("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/gemcitabineVsNaive_scRNAseq.csv")
gemcitabineVsNaive <- gemcitabineVsNaive %>% filter(p_val_adj < 0.001, abs(avg_log2FC) >= 2.5)
gemcitabineVsNaive$gene <- rownames(gemcitabineVsNaive)
write_csv(gemcitabineVsNaive, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/gemcitabineVsNaive_scRNAseq.csv")

sotarasibVsGemcitabine <- FindMarkers(seuratObject, ident.1 = "Sot")
sotarasibVsGemcitabine <- read_csv("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/sotarasibVsGemcitabine_scRNAseq.csv")
sotarasibVsGemcitabine <- sotarasibVsGemcitabine %>% filter(p_val_adj < 0.001, abs(avg_log2FC) >= 2.5)
sotarasibVsGemcitabine$gene <- rownames(sotarasibVsGemcitabine)
write_csv(sotarasibVsGemcitabine, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/sotarasibVsGemcitabine_scRNAseq.csv")

#Now comparing the bulkseq DE genes with scRNAseq DE genes
deSeqBulkRnaSeq <- read_csv("subiaHanxiaoDataAnalysis/extractedData/deSeqBulkRnaSeq_withoutOneGem.csv")
colnames(deSeqBulkRnaSeq) <- c("ensemblID", "baseMean", "LFC-Bulk", "lfcSE", "stat", "pvalue-bulk", "padj-Bulk", "comparison", "randIndex", "gene_name")
deSeqBulkRnaSeq <- deSeqBulkRnaSeq %>% dplyr::select(c("LFC-Bulk","padj-Bulk","comparison","gene_name"))

#First naive vs Sotarasib
sotarasibVsNaiveBulk <- deSeqBulkRnaSeq %>% filter(comparison == "Naive_vs_Sotorasib", `padj-Bulk` < 0.001)
commonGenessotarasibVsNaive <- merge(sotarasibVsNaiveBulk, sotarasibVsNaive, by.x = "gene_name", by.y = "gene")
pValueSotarasib <- randomTesting(Features(seuratObject), sotarasibVsNaive$gene, sotarasibVsNaiveBulk$gene_name)
write_csv(commonGenessotarasibVsNaive, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/commonGenesSotarasibVsNaive_New.csv")

#Next naive vs gemcitibine
GemcitabineVsNaiveBulk <- deSeqBulkRnaSeq %>% filter(comparison == "Naive_vs_Gemcitabine", `padj-Bulk` < 0.001)
commonGenesGemcitibineVsNaive <- merge(GemcitabineVsNaiveBulk, gemcitabineVsNaive, by.x = "gene_name", by.y = "gene")
pValueGemcitabine <- randomTesting(Features(seuratObject), gemcitabineVsNaive$gene, GemcitabineVsNaiveBulk$gene_name)
write_csv(commonGenesGemcitibineVsNaive, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/commonGenesGemcitibineVsNaive_New.csv")

#Sotarasib vs gemcitabine
sotarasibVsGemcitabineBulk <- deSeqBulkRnaSeq %>% filter(comparison == "Sotorasib_vs_Gemcitabine", `padj-Bulk` < 0.001)
commonGenesSotarasibVsGemcitabine <- merge(sotarasibVsGemcitabineBulk, sotarasibVsGemcitabine, by.x = "gene_name", by.y = "gene")
pValueSotarasibGemcitabine <- randomTesting(Features(seuratObject), sotarasibVsGemcitabine$gene, sotarasibVsGemcitabineBulk$gene_name)
write_csv(commonGenesSotarasibVsGemcitabine, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/commonDeGenesBoth/commonGenesSotarasibVsGemcitabine_New.csv")

library(ggvenn)
data <- list(`bulk RNAseq` = GemcitabineVsNaiveBulk$gene_name,
             `scRNAseq` = gemcitabineVsNaive$gene )
venn <- ggvenn(data, c("bulk RNAseq", "scRNAseq"), show_percentage = FALSE, fill_color = c("orange", "pink"))

# Add title using labs
venn <- venn + labs(title = "DE genes between Naive vs Gemcitabine") +
  theme(plot.title = element_text(hjust = 0.5))  
venn
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/vennDiagrams/gemcitabineVsNaive_New.png", venn)


sotarasib <- list(`bulk RNAseq` = sotarasibVsNaiveBulk$gene_name,
                         `scRNAseq` = sotarasibVsNaive$gene )
venn <- ggvenn(sotarasib, c("bulk RNAseq", "scRNAseq"), show_percentage = FALSE, fill_color = c("orange", "pink"))

# Add title using labs
venn <- venn + labs(title = "DE genes between Naive vs Sotarasib") +
  theme(plot.title = element_text(hjust = 0.5))  
venn
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/vennDiagrams/sotarasibVsNaive_New.png", venn)

sotarasibVsGemci <- list(`bulk RNAseq` = sotarasibVsGemcitabineBulk$gene_name,
             `scRNAseq` = sotarasibVsGemcitabine$gene )
venn <- ggvenn(sotarasibVsGemci, c("bulk RNAseq", "scRNAseq"), show_percentage = FALSE, fill_color = c("orange", "pink"))

# Add title using labs
venn <- venn + labs(title = "DE genes between Sotarasib vs Gemcitabine") +
  theme(plot.title = element_text(hjust = 0.5))  
venn
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/vennDiagrams/sotarasibVsGemcitabine_New.png", venn)


########################################################################################
#Barcoded cells
#Combining RNA-sq data with barcode data
singletCells <- read_csv("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/singletList.csv")

singletCells <- singletCells %>%
  mutate(across(everything(), ~ gsub("sotA", "sot", .x)))

colnames(singletCells) <- c("cellID", "barcode") 
collapsedData <- aggregate(barcode ~ cellID, data = singletCells, FUN = function(x) paste(x, collapse = ", "))

seuratCells <- data.frame(cellID = Cells(seuratObject))
commonCellIDList <- merge(seuratCells, collapsedData, by = "cellID")

barcodeSeuratObject <- subset(seuratObject, cells = commonCellIDList$cellID)
subsetSingletCells <- collapsedData[collapsedData$cellID %in% commonCellIDList$cellID, ]
rownames(subsetSingletCells) <- subsetSingletCells$cellID
barcodeSeuratObject$barcode <- subsetSingletCells$barcode

saveRDS(barcodeSeuratObject, "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/rdsObjects/barcodeSeuratObjectHarmony.rds")

#Plot UMAPs

