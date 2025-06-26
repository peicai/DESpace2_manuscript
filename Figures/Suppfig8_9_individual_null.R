rm(list = ls())
library(tidyverse)
path_fig <- "./figure/"
source("./figure_theme.R", echo=TRUE)
revsqrt_trans <- scales::trans_new(
  name = "revsqrt",
  transform = function(x) -sqrt(1 - x),
  inverse = function(x) 1 - x^2
)
names(colors_method) <- c(
  "DESpace_domain","DESpace2",
  "FindMarkers_domain", "FindMarkers (Seurat)", 
  "findMarkers_domain", "findMarkers (scran)",
  "pseudo-bulk_FindMarkers_domain","pseudo-bulk FindMarkers",
  "spatialLIBD_domain", "spatialLIBD")
colors_method["theoretical minPval"] <- "black"

process_and_plot <- function(path_data, dataset_name) {
  # Load Data
  res <- read.csv(paste0(path_data,"results.csv"), row.names = 1)
  if(dataset_name == "ARTISTA"){
    res <- res %>% filter(method != "DESpace_GT")
    if (!grepl("_Banksy$", res$method)) {
      res$method <- paste0(res$method, "_Banksy")
    }
  }
  res <- res %>%
    filter(scenario %in% c("NULL1", "NULL2")) %>%
    separate(method, into = c("method", "annotation"), sep = "_(?=[^_]+$)") %>%
    mutate(
      method = case_when(
        method == "pseudo.bulk_FindMarkers" ~ "pseudo-bulk FindMarkers",
        method == "Theoretical_Min_Pvalue" ~ "theoretical minPval",
        method == "DESpace" ~ "DESpace2",
        method == "FindMarkers" ~ "FindMarkers (Seurat)",
        method == "findMarkers" ~ "findMarkers (scran)",
        TRUE ~ method
      )
    ) %>% filter(annotation %in% c("BayesSpace", "Banksy"))
  p <- ggplot(data = res, aes(x = pval, y = ..ndensity.., lty = domain,
                                              col = method, fill = method)) +
      geom_density(adjust = 0.5, size = 0.3, alpha = 0.1) +
      facet_wrap(scenario~ method, ncol = 5) +
      scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
      guides(col = FALSE,
             lty = guide_legend(ncol = 5, order = 2),
             fill = FALSE) +
      scale_x_continuous("p-value", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
      scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
      .prettify("bw") + 
      my_theme +
      theme(legend.position="bottom",
            legend.direction = "horizontal",
            legend.box="vertical",
            aspect.ratio = 0.65,
            strip.text = element_text(size = 3),
            legend.text=element_text(size=6),
            legend.title=element_text(size=6),
            axis.text.x = element_text(size = 5),
            axis.text.y=element_text(size=5),
            axis.title.y = element_text(size=6),
            axis.title.x = element_text(size=6)) +
      scale_colour_manual(values = colors_method) +
      scale_fill_manual(values = colors_method) 
  rm(res)
  return(p)
}

# Apply Function to Both Datasets
AA1 <- process_and_plot("./data/LIBD_individual/", "LIBD")
AA2 <- process_and_plot("./data/ARTISTA_individual/", "ARTISTA")

ggsave(filename = paste0(path_fig, 'supp_individual_null_LIBD.pdf'),
       plot = AA1,
       device = "pdf",
       width = 7.5,
       height = 3,
       units = "in",
       dpi = 300,
       limitsize = TRUE)
ggsave(filename = paste0(path_fig, 'supp_individual_null_ARTISTA.pdf'),
       plot = AA2,
       device = "pdf",
       width = 7.5,
       height = 3,
       units = "in",
       dpi = 300,
       limitsize = TRUE)