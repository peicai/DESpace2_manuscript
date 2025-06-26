source("figure_theme.R", echo=TRUE)
spe = readRDS("./highly_abundant/ARTISTA_main/combined_simulated_updated.rda")

group1_samples <- c("2DPI_1","5DPI_2", "10DPI_1", "20DPI_2")
group2_samples <- c("2DPI_2","5DPI_1", "10DPI_2", "20DPI_1" )

scenarios <- c("DSP1", "DSP2", "NULL1", "NULL2")
cluster_lists <- list(
  c("0", "NULL"),
  c("0", "2"),
  c("NULL", "NULL"),
  c("0", "0")
)
spe$sample_id <- spe$Sample
sub_sce <- spe
sample_ids <- unique(spe$sample_id)
sorted_ids <- sample_ids[order(as.numeric(sub("DPI.*", "", sample_ids)), 
                               as.numeric(sub(".*_", "", sample_ids)))]

sub_sce$sample_id <- factor(sub_sce$sample_id, levels = sorted_ids)
sce_all <- lapply(levels(sub_sce$sample_id), function(sample){
  sce <- sub_sce[, sub_sce$sample_id == sample]
  sce <- logNormCounts(sce)
  gc()
  sce
})
lapply(sce_all, function(sce) print(unique(sce$sample_id)))
reorder_sce_all <- sce_all[c(1,4,5,8,2,3,6,7)]
output_dir <- "./figure/"
plots_all <- lapply(seq_along(scenarios), function(i){
  scenaio = scenarios[i]
  cluster_list = cluster_lists[[i]]
  plots1 <- lapply(reorder_sce_all, function(spe_sub){
    clust <- ifelse(unique(spe_sub$sample_id) %in% group1_samples,
                    cluster_list[1],
                    cluster_list[2])
    df <- colData(spe_sub) %>% as.data.frame
    df$layer <- .cluster_label(spe = spe_sub,
                               cluster_list = clust,
                               cluster_col = "Cluster")
    
    p <- ggplot(df, aes(x = imagerow, y = imagecol, 
                        colour = layer)) +
      geom_point(size = 0.5) +
      scale_colour_manual(values = c("in" = "darkblue",
                                     "out" = "lightgrey"
      ),
      labels = c("highly abundant", "lowly abundant"),              
      guide = guide_legend(override.aes = list(size = 3))
      )   +
      theme_classic() +
      theme(
        legend.position = "right",
        axis.text = element_blank(),   
        axis.ticks = element_blank(),
        axis.title = element_blank(),  
        axis.line = element_blank()  
      )+
      labs(title = unique(spe_sub$sample_id)) 
    p <- apply_transformations(p, unique(spe_sub$sample_id))
    p
  })
})
flat_plots1 <- unlist(plots_all, recursive = FALSE)
plots <- lapply(reorder_sce_all, function(sce_one) {
  annotation_plot(sce_one,cluster = "Banksy_smooth",
                  xcoord = "sdimx", ycoord = "sdimy")+
    labs(title = unique(sce_one$sample_id)) 
})
flat_plots1 <- c(plots, flat_plots1)

#######################################################
plots_per_page <- 8
total_pages <- 5

fourth_plots <- vector("list", total_pages)
eighth_plots <- vector("list", total_pages)

for (page_i in 1:total_pages) {
  start_index <- (page_i - 1) * plots_per_page
  
  p4 <- flat_plots1[[start_index + 4]] + theme(plot.title = element_blank())
  p8 <- flat_plots1[[start_index + 8]] + theme(plot.title = element_blank())
  
  fourth_plots[[page_i]] <- p4
  eighth_plots[[page_i]] <- p8
}

first_row <- ggpubr::ggarrange(plotlist = fourth_plots, ncol = 5, nrow = 1, legend = "none")
second_row <- ggpubr::ggarrange(plotlist = eighth_plots, ncol = 5, nrow = 1, legend = "none")

combined_rows <- ggpubr::ggarrange(
  first_row,
  second_row,
  ncol = 1,
  nrow = 2,
  heights = c(1, 1),
  legend = "none"
)

# extract the legend
legend <- ggpubr::get_legend(fourth_plots[[2]] + theme(legend.position = "bottom",
                                                       legend.title = element_blank(),
                                                       legend.text = element_text(size = 30)) +
                               guides(
                                 color = guide_legend(override.aes = list(size = 7)))
                             )

# merge plot and legend
final_plot <- ggpubr::ggarrange(
  combined_rows,
  legend,
  ncol = 1,
  heights = c(10, 1)  
)

final_plot <- final_plot + plot_annotation(title = "") &
  theme(plot.title = element_text(size = 20))

ggsave(filename = paste0(output_dir, "Simulation_examples.png"),
       plot = final_plot,
       width = 25,
       height = 10,
       dpi = 300) 
