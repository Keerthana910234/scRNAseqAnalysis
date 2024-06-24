library(ggplot2)
library(dplyr)
library(gridExtra)

# Get unique identifiers
ident_values <- unique(umapDataAllBarcodesAllCount$orig.ident)
ident_values <- "Naive-2"
ident_values <- "Sot"
ident_values <- "Gem"


# Initialize a list to store all combined plots (one combined plot per ident)
combined_plots_list <- list()

for (ident in ident_values) {
  # Get the unique barcodes for this ident
  barcode_values <- unique(umapDataAllBarcodesAllCount$BarcodeNumber[umapDataAllBarcodesAllCount$orig.ident == ident & umapDataAllBarcodesAllCount$sizePerBarcodePerSample > 1])
  
  # Initialize a list to store plots for this ident
  plot_list <- list()
  
  # Iterate over each barcode
  for (barcodeNumber in barcode_values) {
    # Create the plot for each barcode
    barcode_data <- umapDataAllBarcodesAllCount %>%
      filter(orig.ident == ident, BarcodeNumber == barcodeNumber)
    size_sample <- max(barcode_data$sizePerBarcodePerSample, na.rm = TRUE)
    
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
        panel.grid.minor = element_blank()) + facet_wrap(~ orig.ident)  +
      scale_color_identity()  # Use actual color names specified in `aes()`
    
    # Store the plot in the list
    plot_list[[barcodeNumber]] <- p
  }
  
  # Combine the plots for this ident into a single plot using gridExtra
  
  combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 3)) # Adjust ncol as necessary
  combined_plots_list[[ident]] <- combined_plot
}

# Now, combined_plots_list contains all the combined plots for each orig.ident,
# you can access each combined plot via combined_plots_list[[ident_index]] and save or display them as needed

# Optionally, display or save one of the combined plots, for example:
print(combined_plots_list[[1]])  # Replace 1 with appropriate index based on your interest

# Save the combined plot
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/cloneWise/GemCloneWiseNew2.png", plot = combined_plots_list[[1]], width = 30, height = 60, limitsize = F, units = "in")
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/cloneWise/NaiveCloneWiseNew.png", plot = combined_plots_list[[1]], width = 30, height = 75, limitsize = F, units = "in")
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/cloneWise/SotCloneWise.png", plot = combined_plots_list[[1]], width = 30, height = 30, limitsize = F, units = "in")


barcodeNumber = 68
identNew = "Naive-2"
barcodeNew <- umapDataAllBarcodesAllCount %>%
  filter(BarcodeNumber == barcodeNumber) %>%
  filter(orig.ident == identNew)
barcodeNew <- umapDataAllBarcodesAllCount
sizeSample <- max(barcodeNew$sizePerBarcodePerSample, na.rm = TRUE)
shuffled_data <- umapDataAllBarcodesAllCount %>%
  group_by(orig.ident) %>%
  mutate(row_id = row_number()) %>%  # Create an identifier within each group
  ungroup() %>%
  arrange(sample(row_id)) %>%  # Shuffle based on the row identifier
  select(-row_id)
sampleWise <- shuffled_data %>% filter(orig.ident == "Naive-2")
plotSingle <- ggplot(data = shuffled_data, aes(x = umap_1, y = umap_2)) +
  geom_point(color = "lightgrey", alpha = 0.2, size = 2) +
  geom_point(data = sampleWise, aes(x = umap_1, y = umap_2, color = "#00ba37"), size = 2, alpha = 1) +  # Ensure this layer is on top
  # labs(title = paste("UMAP of", identNew, "-", barcodeNumber, "size", sizeSample)) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.6) +
  scale_color_identity()  # Use actual color names specified in `aes()`
plotSingle
ggsave("/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/cloneWisePlotsFinal/Naive.png", height = 8, width = 12)
svglite(filename = "/projects/b1042/GoyalLab/Keerthana/subiaHanxiaoDataAnalysis/plots/cloneWisePlotsFinal/Naive.svg", width = 12, height = 8)
plot(plotSingle)
dev.off()

library(svglite)
