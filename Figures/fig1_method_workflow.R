rm(list = ls())
library(scuttle)
library(ggplot2)
library(SpatialExperiment)
library(dplyr)
path <- "./"
load(paste0(path, "Results/DESpace_Banksy.rda"))
load(paste0(path, "Results/DESpace_GT.rda"))
load(paste0(path, "Banksy_harmony_smooth_rotated.rda"))
gc()
## Visualization
##################### plotting scheme #####################
colData(spe) <- cbind(colData(spe), spatialCoords(spe))
spe <- logNormCounts(spe)
sample_ids <- unique(spe$sample_id)
# rotate tissue coordinates to ensure all sections look consistent.
rotate_coord <- function(spe, sample) {
  spe_object <- spe[, spe$sample_id == sample]
  sdimx <- spe_object[["sdimx"]]
  sdimy <- spe_object[["sdimy"]]
  if (sample %in% c('2DPI_1', '5DPI_3')) {
    # Rotate 180 degrees: Reverse both x and y
    sdimx <- -sdimx
    sdimy <- -sdimy
    # Rotate 90 degrees counterclockwise: Swap x and y
    sdimx <- -sdimy
    sdimy <- -spe_object[["sdimx"]]
  } else if (sample %in% c('5DPI_1', '5DPI_2', '10DPI_2', 
                           '15DPI_1', '15DPI_2', '15DPI_3',
                           '15DPI_4', '20DPI_1', '20DPI_2')) {
    # # Rotate 90 degrees counterclockwise: Swap x and y
    # sdimx <- -sdimy
    # sdimy <- spe_object[["sdimx"]]
  } else if (sample %in% c('10DPI_3')) {
    # Rotate 90 degrees clockwise: Swap x and y, reverse new y
    sdimx <- sdimy
    sdimy <- -spe_object[["sdimx"]]
  } else if (sample %in% c('2DPI_2', '10DPI_1')) {
    sdimx <- -sdimy
    sdimy <- spe_object[["sdimx"]]
  }
  
  # Scale to ensure all coordinates are positive
  spe_object$sdimx <- sdimx - min(sdimx) + 1
  spe_object$sdimy <- sdimy - min(sdimy) + 1
  spe_object$sdimx <- -spe_object$sdimx
  return(spe_object)
}
all_spe <- lapply(sample_ids, function(x) rotate_coord(spe, x))
spe <- do.call(cbind, all_spe); rm(all_spe)
# save(spe, file = paste0(path, "Banksy_harmony_smooth_rotated.rda"))
# Define a function to create plots for each gene
create_plot <- function(sce_one, gene, gene_title) {
  sce_one$sample_id <- factor(sce_one$sample_id, levels = c("2DPI_1", "2DPI_2",
                                            "5DPI_1", "5DPI_2",
                                            "10DPI_1", "10DPI_2",
                                            "20DPI_1", "20DPI_2"))
  df <- data.frame(Exp. = logcounts(sce_one)[gene, ],
                   row = sce_one$sdimx,
                   col = sce_one$sdimy,
                   sample = sce_one$sample_id)
  
  plot_title <- gene_title
  filter <- c()
  (p <- ggplot(df, aes(x = row, y = col, colour = Exp.)) + 
      geom_point(alpha = 0.9, size = 0.1 ) +
      scale_colour_gradientn(colors = c("lightgrey", "darkblue"),
                             breaks = c(0, max(df$Exp.)),  # Dynamic breaks
                             labels = c("low", "high")) +
      theme_bw() +
      scale_y_continuous(expand = c(0, 0)) +
      facet_wrap(~sample, scales = "free", nrow = 2) +
      theme(legend.position = "left",
            title = element_text(colour = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(color = "white",size = 15),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color  =  NA),
            legend.background = element_rect(fill = "white"),
            legend.text = element_text(color = "white", hjust = 0,size = 12),
            legend.justification = c("center", "top"),
            legend.direction = "vertical",
            legend.title = element_text(color = "white",hjust = 0,size = 15),  # Center the title
            plot.background = element_rect(color = "white", fill = "white")
      ) + guides(color = guide_colourbar(theme = theme(
        legend.text.position = "left"
      ))) 
    # +
    #labs(title = plot_title) 
  )
}
######################### plot ###########################
load(paste0(path, "RegionsGT_DSP1_simulated.rda"))
gene_id <- names(inf)[sample(length(inf), 25)]
plots1 <- lapply(gene_id, function(each_gene){
  create_plot(spe, each_gene, gene_title=each_gene)
})
flat_plots1 <- unlist(plots1, recursive = FALSE)

output_dir <- "./Results/SVG/"
method_name = "DESpace"

library(patchwork)
pdf(paste0(output_dir, "DSP2_GeneExpr_lowest5Pval_", method_name, ".pdf"), width = 36, height = 9)
combined_plot1 = ggpubr::ggarrange(plotlist = plots1, ncol=2)
print(combined_plot1)
dev.off()

pdf(paste0(output_dir, "DSP2_GeneExpr_lowest1Pval_", method_name, ".pdf"), width = 9, height = 4.5)
print(plots1[[1]])
dev.off()

# generate a list of FeaturePlots
feature = gene_id[1]
plots <- lapply(c("20DPI_1","5DPI_1",
                  "20DPI_2","5DPI_2"), function(sample_id) {
  # Subset spe for each sample
  spe_j <- spe[, colData(spe)$sample_id == sample_id]
  # Create FeaturePlot for the sample
  plot <- DESpace::FeaturePlot(spe_j, feature, 
                      cluster_col = "BANKSY_snn_res0.3_smooth",
                      coordinates = c("sdimx", "sdimy"), cluster = 'Layer1',
                      platform = "Stereo-seq",
                      diverging = TRUE,
                      #sf_dim = 50,
                      point_size = 0.2,
                      low = "lightgrey",
                      high = "darkblue",
                      linecolor = "#005682",
                      linewidth = 0.6) +
    theme(legend.position = "right",
          legend.key.size = unit(0.5, 'cm')) +
    labs(color = "") #+ ggtitle(sample_id) 
  
  return(plot)
})
combined_plot <- wrap_plots(plots, ncol = 2) + 
  # common legend
  plot_layout(guides = 'collect')  

pdf(paste0(output_dir, "DSP2_GeneExprInd_lowest1Pval_", method_name, ".pdf"), width = 8, height = 6)
combined_plot
dev.off()