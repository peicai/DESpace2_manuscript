rm(list=ls())
library(dplyr)
library(ggplot2)

################################ ARTISTA ###############################
library(muSpaData)
spe <- Wei22_full()
df <- as.data.frame(colData(spe))
p <- df %>% 
  ggplot(aes(x = sdimx, y = sdimy, colour = Banksy_smooth)) +
  geom_point(size = 0.1) +
  facet_wrap(~ sample_id) +
  theme_classic() + 
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),    
    axis.ticks = element_blank(),   
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"))+
  guides(colour = guide_legend(override.aes = list(size = 2)))

path_fig <- "./figure/"
ggsave(filename = paste0(path_fig, 'supp_ARTISTA_Banksy.pdf'),
       plot = p,
       device = "pdf",
       width = 11,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
rm(spe)
############################### LIBD ######################################
library(SingleCellExperiment)
library(patchwork)
data_path <- "./LIBD_filtered/"
sample_id <- list(c("151507", "151509"),
                  c("151669", "151671"),
                  c("151673", "151675"))

sce_objects <- lapply(unlist(sample_id), function(id) {
  print(id)
  load(paste0(data_path, id, "_sce_LIBD.rda"))
  sce_one$layer_guess_reordered_droplevel <- sce_one$layer_guess_reordered
  ind <- sce_one$layer_guess_reordered_droplevel %in% c("Layer1", "Layer2")
  if (sum(ind) != 0) sce_one[,ind]$layer_guess_reordered_droplevel <- "Layer6"
  sce_one$layer_guess_reordered_droplevel <- droplevels(sce_one$layer_guess_reordered_droplevel)
  sce_one <- sce_one[,!is.na(sce_one$layer_guess_reordered_droplevel)]
  print(levels(sce_one$layer_guess_reordered_droplevel))
  return(sce_one)
})
common_gene <- Reduce(intersect, lapply(sce_objects, rownames))
sce_objects_common <- lapply(sce_objects, function(sce) {
  colnames <- c("sample_name", "imagerow", "imagecol", "col", "row",
                "layer_guess_reordered_droplevel",
                "layer_guess_reordered")
  colData(sce) <- colData(sce)[, colnames]
  colnames(colData(sce)) <- c("sample_name", "pxl_row_in_fullres", 
                     "pxl_col_in_fullres", "array_col", "array_row",
                     "layer", "original")
  sce[common_gene, , drop = FALSE]
})

custom_colors <- c("#FF7F00", "black", "darkgrey", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
all_clusters <- c("WM","Layer1", "Layer2", "Layer3", "Layer4", "Layer5", "Layer6")  

plots <- lapply(seq_along(sce_objects_common), function(i) {
  sce <- sce_objects_common[[i]]
  sce$original <- factor(sce$original, levels = all_clusters)
  
  p <- BayesSpace::clusterPlot(sce, label = "original", 
                               palette = NULL, size = 0.05, color=NA) +
    scale_fill_manual(
      values = custom_colors,
      name = "",  
      labels = all_clusters 
    ) +
    labs(title = unique(sce$sample_name)) +
    guides(fill = guide_legend(ncol = 4)) +
    theme(plot.title = element_text(hjust = 0.5, size = 4),
          legend.text = element_text(size = 3.5),
          legend.key.size = unit(0.4, "cm"))
  
  # Hide legend for all but first plot
  if (i != 1) p <- p + guides(fill = "none")
  p
})

(combined_plot <- wrap_plots(plots, ncol = 3) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")) 

plots2 <- lapply(sce_objects_common, function(sce) {
     sce$layer <- factor(sce$layer, levels = c("WM","Layer3", "Layer4", "Layer5", "Layer6"))
     BayesSpace::clusterPlot(sce, label = "layer", palette = NULL, size = 0.05, color=NA) +
         scale_fill_manual(
             values = custom_colors[c(1,4:7)],
             name = "",  
             labels = c("WM","Layer3", "Layer4", "Layer5", "Layer6 (merged L1 and 2)")  
           ) +
         labs(title = unique(sce$sample_name)) + 
    guides(fill = guide_legend(ncol = 3)) +
       theme(plot.title = element_text(hjust = 0.5, size = 4),
             legend.text = element_text(size = 3.5),
             legend.key.size = unit(0.4, "cm"))
   })
(combined_plot2 <- wrap_plots(plots2, ncol = 3) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom"))

p <- cowplot::plot_grid(combined_plot, combined_plot2, 
                        ncol = 1, rel_widths = c(1,1))
path_fig <- "./figure/"
ggsave(filename = paste0(path_fig, 'supp_LIBD_annotations.pdf'),
       plot = p,
       device = "pdf",
       width = 3.5,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)