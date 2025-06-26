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
my_theme = theme(legend.position = "none", aspect.ratio = 0.8,
                 plot.title = element_text(hjust = 0, face = "bold", size=13),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=12),
                 axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),  
                 axis.text.y = element_text(size = 11),
                 axis.title.x = element_text(size = 11) ,
                 axis.title.y = element_text(size = 11)
                 
)
process_data <- function(gene_list, scenario_name) {
  sub_pval <- pval[rownames(pval) %in% unlist(gene_list), ]
  
  long_df <- sub_pval %>%
    rownames_to_column("gene") %>%
    pivot_longer(cols = -gene, names_to = "Method", values_to = "pval") %>%
    mutate(scenario = scenario_name)
  
  long_df %>%
    separate(Method, into = c("method", "annotation"), sep = "_(?=[^_]+$)")
}
generate_uniform_res <- function(seed, scenario, cluster_type = "BayesSpace", method_name = "Theoretical_Min_Pvalue") {
  set.seed(seed)  # Set seed for reproducibility
  uniform_pvals <- matrix(runif(6 * 1e6, min = 0, max = 1), nrow = 1e6, ncol = 6)
  uniform_min_pvals <- apply(uniform_pvals, 1, min)
  n <- length(uniform_min_pvals)
  
  # Create the result data frame
  res_uniform <- data.frame(
    pval = uniform_min_pvals,
    #padj = uniform_min_pvals,
    method = rep(method_name, n),
    annotation = rep(cluster_type, n),
    gene = as.character(1:n),
    scenario = rep(scenario, n)
  )
  
  return(res_uniform)
}

process_and_plot <- function(path_data, dataset_name) {
  # Load Data
  pval <- read.csv(paste0(path_data, "pval_all.csv"), row.names = 1)
  gene_lists <- readRDS(paste0(path_data, "geneGT.rds"))
  
  # Process Data
  res1 <- process_data(gene_lists[3], "NULL1")
  res2 <- process_data(gene_lists[4], "NULL2")
  res_uniform1 <- generate_uniform_res(seed = 123, scenario = "NULL1")
  res_uniform2 <- generate_uniform_res(seed = 234, scenario = "NULL2")
  
  res <- rbind(res1, res2, res_uniform1, res_uniform2) %>%
    mutate(
           method = case_when(
             method == "pseudo.bulk_FindMarkers" ~ "pseudo-bulk FindMarkers",
             method == "Theoretical_Min_Pvalue" ~ "theoretical minPval",
             method == "DESpace" ~ "DESpace2",
             method == "FindMarkers" ~ "FindMarkers (Seurat)",
             method == "findMarkers" ~ "findMarkers (scran)",
             TRUE ~ method
           ),
           # annotation = case_when(
           #   annotation %in% c("BayesSpace", "Banksy") ~ "domain",
           #   annotation == "GT" ~ "GT",
           #   TRUE ~ method
           # ),
           allmethod = paste0(method, "_", annotation)
    ) %>% filter(annotation %in% c("BayesSpace", "Banksy"))
  # Generate Plot
  p <- ggplot(res, aes(x = pval, y = ..ndensity.., col = method, fill = method)) +
    geom_density(adjust = 0.5, size = 0.3, alpha = 0.5) +
    facet_grid(scenario ~ method) +
    scale_alpha_manual(values = c("TRUE" = 0.1, "FALSE" = 0.4)) +
    scale_x_continuous("FDR", breaks = seq(0, 1, 0.2), expand = c(0, 0.05)) +
    scale_y_continuous("density", breaks = c(0, 1), expand = c(0, 0.1)) +
    .prettify("bw") + 
    my_theme +
    scale_colour_manual(values = colors_method) +
    scale_fill_manual(values = colors_method) +
    guides(color = "none", fill = "none") +
    theme(legend.position = "none", strip.text = element_text(size = rel(1.2))) +
    labs(title = dataset_name)
  rm(pval, res, gene_lists, res1, res2, res_uniform1, res_uniform2)
  return(p)
}

# Apply Function to Both Datasets
AA1 <- process_and_plot("./data/LIBD_main/", "LIBD")
AA2 <- process_and_plot("./data/ARTISTA_main/", "ARTISTA")

# Combine Plots
(AA <- ggpubr::ggarrange(AA1, AA2, heights = c(1, 1), ncol = 1))

ggsave(filename = paste0(path_fig, 'supp_global_null.pdf'),
       plot = AA,
       device = "pdf",
       width = 11,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)