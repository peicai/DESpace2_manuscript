rm(list = ls())
source("figure_theme.R", echo=TRUE)
library(DESpace) # require R >= 4.5
library(patchwork)
library(cowplot)
spe = readRDS("./Data/highly_abundant/ARTISTA_main/simulated_sub.rda")

group1_samples <- c("2DPI_1","5DPI_2", "10DPI_1", "20DPI_2")
group2_samples <- c("2DPI_2","5DPI_1", "10DPI_2", "20DPI_1" )

features <- c("AMEX60DD047626", "AMEX60DD026790", "AMEX60DD020041", "AMEX60DDU001008168")
scenarios <- c("DSP1", "DSP2", "NULL1", "NULL2")

pp <- list()

for (i in seq_along(features)) {
  feature <- features[i]
  scenario_title <- scenarios[i]
  
  plots <- lapply(spe, function(spe_j) {
    p <- FeaturePlot(spe_j, feature,
                              coordinates = c("sdimx", "sdimy"),
                              platform = "Stereo-seq", ncol = 1,
                              diverging = TRUE,
                              low = "grey85", mid = "firebrick", high = "red",
                              point_size = 0.1, legend_exprs = TRUE) + 
      theme(legend.position = "none",
            legend.key.size = unit(0.5, 'cm')) +
      labs(color = "") + 
      ggtitle(unique(spe_j$sample_id)) +
      theme(aspect.ratio = 1,
            plot.title = element_text(size = 10))
    return(p)
  })
  combined_plot <- wrap_plots(plots, ncol = 8) 
  pp[[i]] <- combined_plot
}

# Add title to each combined plot individually
for (i in seq_along(pp)) {
  pp[[i]] <- pp[[i]] + plot_annotation(title = scenarios[i]) &
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
}

# Stack vertically with cowplot
final_plot <- plot_grid(plotlist = pp, ncol = 1, align = 'v', axis = 'lr')

ggsave("supp_Simulation_Expression.pdf", plot = final_plot, width = 25, height = 15, units = "in", dpi = 300)