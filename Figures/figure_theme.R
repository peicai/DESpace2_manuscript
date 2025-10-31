library(iCOBRA)
library(RColorBrewer)
library(ggplot2)
library(SingleCellExperiment)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(pROC)
library(cowplot)
library(purrr) 
library(scales)
library(scuttle)
library(patchwork)
library(tidyr)

all_colours = c(
  "#08306B",
  "#A6CEE3",#"#1F78B4", 
  
  "#33A02C",
  "#B2DF8A",
  
  "gold",
  "#F0E68C",
  
  "plum2",
  "mediumpurple",
  
  "#FFA07A",
  "#B22222",
  "white"
)
names(all_colours) <- c(
  "DESpace_GT", "DESpace_BayesSpace",
  "FindMarkers_GT", "FindMarkers_BayesSpace",
  "findMarkers_GT", "findMarkers_BayesSpace",
  "pseudo.bulk_FindMarkers_GT", "pseudo.bulk_FindMarkers_BayesSpace",
  "spatialLIBD_GT","spatialLIBD_BayesSpace")
colors_method <- all_colours
names(colors_method) <- c(
  names(all_colours)[1:6],
  "pseudo-bulk_FindMarkers_GT","pseudo-bulk_FindMarkers_BayesSpace",
  names(all_colours)[9:10])
.prettify <- function(theme = NULL, ...) {
  if (is.null(theme)) theme <- "classic"
  base <- paste0("theme_", theme)
  base <- getFromNamespace(base, "ggplot2")
  base(base_size = 8) + theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, color = "lightgrey"),
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text = element_text(color = "black"),
    legend.key.size = unit(2, "mm"),
    strip.background = element_rect(fill = NA),
    plot.margin = unit(rep(1, 4), "mm"),
    panel.spacing = unit(0, "mm"),
    legend.margin = margin(0,0,1,0,"mm"),
    ...)}

my_theme <- theme(aspect.ratio = 0.67,
                  #legend.box.just = "left",
                  panel.grid = element_blank(),
                  panel.spacing = unit(1, "mm"),
                  panel.border = element_rect(color = "grey"),
                  strip.text = element_text(size = 4),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(1.5)),
                  axis.text.y=element_text(size=rel(1.5)),
                  axis.title.y = element_text(size=rel(1.5)),
                  axis.title.x = element_text(size=rel(1.5)),
                  axis.title.x.top = element_text(size=rel(1.5)),
                  legend.title=element_text(size=rel(1.5)),
                  legend.text=element_text(size=rel(1.5)),
                  legend.key.width=unit(0.4, "cm"),
                  legend.position="bottom",
                  legend.direction = "horizontal",
                  legend.box="horizontal",
                  strip.text.x = element_text(size = rel(2)),
                  strip.text.y = element_text(size = rel(2)),
                  legend.margin=margin())
output_dir <- "./manuscript_code/figure/"

.cluster_label <-
  function(spe, 
           cluster_list = "NULL", 
           cluster_col = "layer_guess_reordered"){
    metadata <- as.data.frame(colData(spe))
    if(cluster_list == "NULL"){
      one_layer <- sample(c("in", "out"), size = nrow(metadata), replace = TRUE)
      print(table(one_layer))
    }else{
      layer_rename <- c()
      layer_rename <- as.character(metadata[[cluster_col]])
      layer_rename[layer_rename != cluster_list] <- 'Other'
      layer_rename <- as.factor(layer_rename)
      layer_rename <- ifelse(layer_rename == "Other", "out", "in")
      one_layer <- layer_rename
    }
    return(one_layer)
  }

apply_transformations <- function(p, sample_id) {
  if (sample_id %in% c('2DPI_1', '5DPI_3')) {
    p <- p + coord_flip() + scale_x_reverse(expand = c(0, 0)) + scale_y_reverse(expand = c(0, 0))
  } else if (sample_id %in% c('2DPI_2', '2DPI_3', '20DPI_3', '10DPI_1')) {
    p <- p + coord_flip()
  } else if (sample_id %in% c('5DPI_1', '5DPI_2', '10DPI_2', '15DPI_1', '15DPI_2', '15DPI_3', 
                              '15DPI_4', '20DPI_1', '20DPI_2')) {
    p <- p + scale_x_reverse(expand = c(0, 0))
  } else if (sample_id %in% c('10DPI_3')) {
    p <- p + scale_y_reverse(expand = c(0, 0))
  }
  return(p)
}

annotation_plot <- function(sce_one,cluster = "BANKSY_snn_res0.1_smooth",
                            xcoord = "sdimx", ycoord = "sdimy") {
  df <- as.data.frame(colData(sce_one))
  df$layer_guess_merged <- as.character(df[[cluster]])
  df <- df[!is.na(df$layer_guess_merged),  ]
  print(head(df))
  (p <- ggplot(df, aes(x = df[[xcoord]], y = df[[ycoord]], 
                       colour = layer_guess_merged)) + 
      geom_point(size = 0.5 ) +
      theme_void() +  
      scale_color_manual(values = c(
                         # "blue", "black","#1B9E77", "orange",  "firebrick"
                         "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
      ) +
      labs(color = "") +
      guides(colour = guide_legend(override.aes = list(size = 3))) +  # Set legend point size here
      theme(legend.position = "left",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.title = element_text(size = 15),
            legend.text = element_text(color = "white", hjust = 0,size = 12),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size=1)
      )  )
  p <- apply_transformations(p, unique(sce_one$sample_id))
  print(p)
}

revsqrt_trans <- scales::trans_new(
  name = "revsqrt",
  transform = function(x) -sqrt(1 - x),
  inverse = function(x) 1 - x^2
)