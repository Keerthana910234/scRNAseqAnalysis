library(Seurat)
library(ggplot2)
library(dplyr)

#Plotting neat UMAPs
plotUmap <- function(dataframe, xAxis, yAxis, colourBy, legend = "right"){
  # Convert string inputs into symbols and then into expressions
  p <- ggplot(data = dataframe, aes(x = !!rlang::sym(xAxis), 
                                    y = !!rlang::sym(yAxis), 
                                    color = !!rlang::sym(colourBy))) +
    geom_point(alpha = 0.8, size = 1) +  # Adding point layer with some transparency
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

#Plotting neat PCA plots
plotPca <- function(dataframe, xAxis, yAxis, colourBy, legend = "right"){
  # Convert string inputs into symbols and then into expressions
  p <- ggplot(data = dataframe, aes(x = !!rlang::sym(xAxis), 
                                    y = !!rlang::sym(yAxis), 
                                    color = !!rlang::sym(colourBy))) +
    geom_point(alpha = 0.8, size = 1) +  # Adding point layer with some transparency
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

seuratObject <- readRDS("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/extractedData/rdsObjects/barcodeSeuratObjectHarmony.rds")

seuratObject$barcode

umapDataAllBarcodes <- FetchData(seuratObject, vars = c("umap_1", "umap_2", "barcode", "orig.ident"))

umapDataAllBarcodes$cellID <- rownames(umapDataAllBarcodes)


umapDataAllBarcodesCount <- umapDataAllBarcodes %>%
  group_by(barcode, orig.ident) %>%  # Group by barcode and sample
  mutate(sizePerBarcodePerSample = n()) %>%  # Count number of entries per group
  ungroup()  

umapDataAllBarcodesCount$BarcodeNumber <- as.numeric(factor(umapDataAllBarcodesCount$barcode))



barcodeNumber = 68
identNew = "Sot"
barcodeNew <- umapDataAllBarcodesCount %>%
  filter(BarcodeNumber == barcodeNumber) %>%
  filter(orig.ident == identNew)

barcodeNew <- umapDataAllBarcodesAllCount
sizeSample <- max(barcodeNew$sizePerBarcodePerSample, na.rm = TRUE)

barcode_data <- umapDataAllBarcodesAllCount %>%
  filter(orig.ident == ident, BarcodeNumber == barcodeNumber)
size_sample <- max(barcode_data$sizePerBarcodePerSample, na.rm = TRUE)

# Prepare Data

# Define a color mapping for 'orig.ident'
color_mapping <- c("Sot" = "hotpink3", "Gem" = "navy", "Naive-2" = "orange")  # Add mappings as needed


newPlotUmap <- function( barcodeNumber, ident){
  background_data <- umapDataAllBarcodesCount
  highlight_data <- umapDataAllBarcodesCount %>%
    filter(orig.ident == ident, BarcodeNumber == barcodeNumber)
  
  plot <- ggplot() +
    geom_point(data = background_data, aes(x = umap_1, y = umap_2), 
               color = "lightgrey", size = 2, alpha = 0.4) +
    geom_point(data = highlight_data, aes(x = umap_1, y = umap_2, color = orig.ident), 
               size = 2, alpha = 1) +
    scale_color_manual(values = color_mapping) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  print(plot)
  return(plot)
}

newPlotUmap_allIdents <- function( barcodeNumber){
  background_data <- umapDataAllBarcodesCount
  highlight_data <- umapDataAllBarcodesCount %>%
    filter(BarcodeNumber == barcodeNumber)
  
  plot <- ggplot() +
    geom_point(data = background_data, aes(x = umap_1, y = umap_2), 
               color = "lightgrey", size = 2, alpha = 0.4) +
    geom_point(data = highlight_data, aes(x = umap_1, y = umap_2, color = orig.ident), 
               size = 2, alpha = 1) +
    scale_color_manual(values = color_mapping) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  print(plot)
  return(plot)
}

path = "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/cloneWisePlots/"
umapBiggestCloneSot <- newPlotUmap(68, "Sot")
umapBiggestCloneSot
library(svglite)
ggsave(filename = paste0(path, "commonAcrossAll_Sot.png"), umapBiggestCloneSot, height = 8, width = 12)
svglite(filename = paste0(path, "commonAcrossAll_Sot.svg"), width = 12, height = 8)
plot(umapBiggestCloneSot)
dev.off()

umapBiggestCloneGem <- newPlotUmap(68, "Gem")
umapBiggestCloneGem
library(svglite)
ggsave(filename = paste0(path, "commonAcrossAll_Gem.png"), umapBiggestCloneGem, height = 8, width = 12)
svglite(filename = paste0(path, "commonAcrossAll_Gem.svg"), width = 12, height = 8)
plot(umapBiggestCloneGem)
dev.off()

umapBiggestCloneNaive <- newPlotUmap(68, "Naive-2")
umapBiggestCloneNaive
library(svglite)
ggsave(filename = paste0(path, "commonAcrossAll_Naive.png"), umapBiggestCloneNaive, height = 8, width = 12)
svglite(filename = paste0(path, "commonAcrossAll_Naive.svg"), width = 12, height = 8)
plot(umapBiggestCloneNaive)
dev.off()

umapGem1 <- newPlotUmap(567, "Gem")
umapGem1
library(svglite)
ggsave(filename = paste0(path, "umapGem1.png"), umapGem1, height = 8, width = 12)
svglite(filename = paste0(path, "umapGem1.svg"), width = 12, height = 8)
plot(umapGem1)
dev.off()

umapGem2 <- newPlotUmap(143, "Gem")
umapGem2
library(svglite)
ggsave(filename = paste0(path, "umapGem2.png"), umapGem2, height = 8, width = 12)
svglite(filename = paste0(path, "umapGem2.svg"), width = 12, height = 8)
plot(umapGem2)
dev.off()
1281

umapGem3 <- newPlotUmap(1281, "Gem")
umapGem3
library(svglite)
ggsave(filename = paste0(path, "umapGem3.png"), umapGem3, height = 8, width = 12)
svglite(filename = paste0(path, "umapGem3.svg"), width = 12, height = 8)
plot(umapGem3)
dev.off()


umapSot1 <- newPlotUmap(3479, "Sot")
umapSot1
library(svglite)
ggsave(filename = paste0(path, "umapSot1.png"), umapSot1, height = 8, width = 12)
svglite(filename = paste0(path, "umapSot1.svg"), width = 12, height = 8)
plot(umapSot1)
dev.off()

umapSot2 <- newPlotUmap(2138, "Sot")
umapSot2
library(svglite)
ggsave(filename = paste0(path, "umapSot2.png"), umapSot2, height = 8, width = 12)
svglite(filename = paste0(path, "umapSot2.svg"), width = 12, height = 8)
plot(umapSot2)
dev.off()

umapSot3 <- newPlotUmap(3572, "Sot")
umapSot3
library(svglite)
ggsave(filename = paste0(path, "umapSot3.png"), umapSot3, height = 8, width = 12)
svglite(filename = paste0(path, "umapSot3.svg"), width = 12, height = 8)
plot(umapSot3)
dev.off()

umapNaive1 <- newPlotUmap(1519, "Naive-2")
umapNaive1
library(svglite)
ggsave(filename = paste0(path, "umapNaive1.png"), umapNaive1, height = 8, width = 12)
svglite(filename = paste0(path, "umapNaive1.svg"), width = 12, height = 8)
plot(umapNaive1)
dev.off()

umapNaive2 <- newPlotUmap(1732, "Naive-2")
umapNaive2
library(svglite)
ggsave(filename = paste0(path, "umapNaive2.png"), umapNaive2, height = 8, width = 12)
svglite(filename = paste0(path, "umapNaive2.svg"), width = 12, height = 8)
plot(umapNaive2)
dev.off()

umapNaive3 <- newPlotUmap(2076, "Naive-2")
umapNaive3
library(svglite)
ggsave(filename = paste0(path, "umapNaive3.png"), umapNaive3, height = 8, width = 12)
svglite(filename = paste0(path, "umapNaive3.svg"), width = 12, height = 8)
plot(umapNaive3)
dev.off()

umapBiggestClone <- newPlotUmap_allIdents(68)
umapBiggestClone
ggsave(filename = paste0(path, "commonAcrossAll.png"), umapBiggestClone, height = 8, width = 12)
svglite(filename = paste0(path, "commonAcrossAll.svg"), width = 12, height = 8)
plot(umapBiggestClone)
dev.off()






