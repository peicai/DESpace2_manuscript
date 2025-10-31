rm(list = ls())
source("./figure_theme.R", echo=TRUE)
library(dplyr)
library(scales)    
library(colorspace)
path_fig = "./figure/"
path_data <- "./data/"
shape_border = c(0,6,1,2,3)
shape_fill = c(15,25,21,17,3)
sel_shape = c(1,3,2,4,5)
shuffles <- c("0", "001", "005", "01", "02")
fig_theme <- theme(strip.background = element_blank(),
                   strip.text = element_blank(),
                   axis.text.x = element_text(angle = 25, hjust=0.5, size = rel(1.5)),
                   axis.text.y=element_text(size=rel(1.5)),
                   axis.title.y = element_text(size=rel(1.3)),
                   axis.title.x = element_text(size=rel(1.3)),
                   legend.title=element_text(size=rel(1)),
                   legend.text=element_text(size=rel(1)),
                   legend.key.width=unit(1, "cm"),
                   legend.key = element_rect(fill = "grey90", colour = "white"),
                   legend.background = element_blank(),
                   aspect.ratio = 1, legend.position="bottom",
                   legend.box="vertical", legend.margin=margin())
########################################################################
########################## ARTISTA global test ############################
########################################################################
######################### Load data #############################
pval <- read.csv(paste0(path_data,"pval.csv"), row.names = 1)
padj <- read.csv(paste0(path_data,"padj.csv"), row.names = 1)

gene_lists <- readRDS(paste0(path_data, "geneGT.rds"))

scenario <- list(
    gene1 = unlist(gene_lists[c(1:2)]), # DSP1 + DSP2
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  )

names(colors_method) <- gsub("_BayesSpace", "", names(colors_method))
names(all_colours) <- gsub("_BayesSpace", "", names(all_colours))

gg_fdr_ARTISTA <- list()
gene <- c(scenario$gene1, scenario$gene2)
all_colours[["spatialLIBD_pairwise"]] <- all_colours[["spatialLIBD"]]
sub_pval <- pval[rownames(pval) %in% gene, ]
sub_padj <- padj[rownames(padj) %in% gene, ]
truth <- data.frame(status = c(rep(1, length(scenario$gene1)), rep(0, length(scenario$gene2))))
rownames(truth) <- gene
methods_list <- c("DESpace",
                  "FindMarkers", 
                  "findMarkers",
                  "pseudo.bulk_FindMarkers",
                  "spatialLIBD_pairwise")
combo <- expand.grid(method = methods_list, shuffle = shuffles, stringsAsFactors = FALSE)
methods_order <- paste(combo$method, combo$shuffle, sep = "_")
methods_names <- c("DESpace2",
                   "FindMarkers (Seurat)",
                   "findMarkers (scran)",
                   "pseudo-bulk FindMarkers",
                   "spatialLIBD")

grad_list <- lapply(methods_list, function(m){
  base_col <- all_colours[m]        
  pal <- gradient_n_pal(c(lighten(base_col, 0.6), base_col, darken(base_col, 0.6)))
  pal(seq(0, 1, length.out = length(shuffles)))  
})
names(grad_list) <- methods_list

# Flatten to a named vector corresponding to method_shuffle
grad_colors <- unlist(lapply(methods_list, function(m){
  setNames(grad_list[[m]], paste(m, shuffles, sep="_"))
}))

# Create COBRA data and calculate performance
DF_COBRA <- COBRAData(pval = data.frame(sub_pval), 
                      padj = data.frame(sub_padj), 
                      truth = truth)
perf <- calculate_performance(DF_COBRA, binary_truth = "status")

# Prepare for plotting
cobra_plot <- prepare_data_for_plot(
  perf,
  colorscheme = c(grad_colors, "white"),
  incloverall = FALSE,
  facetted = TRUE,
  conditionalfill = FALSE
)

# Plot ROC and FDR/TPR curves
(gg_fdr_ARTISTA_v3 <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                       pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.4), limits = c(0,0.4)) +
    scale_color_manual(values = grad_colors,
                       name = "",
                       breaks=methods_order)+
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(c(rep(0,5), rep(c(1,6),5),
                  rep(2,5), rep(3,5)), 3), 
    stroke = 2, 
    alpha = 1) + 
    scale_fill_manual(values = rep(grad_colors,3),
                      name = "",
                      breaks= methods_order)+
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(c(rep(15,5), rep(c(21,25),5),
                  rep(17,5), rep(3,5)), 3), 
    stroke = 2, alpha = 0.25)+
    fig_theme + theme(legend.position = "none"
    ) )

####################### individual ######################
library(readr)
library(scales)
library(colorspace)
source("./figure_theme.R", echo=TRUE)
result_combined <- read_csv(paste0(path_data,"results.csv"))
scenarios <- list(
  c("DSP1", "NULL1", "NULL2"),
  c("DSP2", "NULL1", "NULL2"),
  c("DSP1", "DSP2", "NULL1", "NULL2")
)
create_wide_df <- function(data, value_col) {
  data %>%
    select(gene, method, !!sym(value_col), domain, scenario) %>% 
    pivot_wider(
      names_from = method,
      values_from = !!sym(value_col)
    ) %>% 
    as.data.frame() %>%
    { 
      rownames(.) <- with(., paste(gene, domain, scenario, sep = "_")) 
      select(., -c(gene, domain, scenario)) 
    }
}

shuffles <- c("0", "001", "005", "01", "02")
q <- 3
shape_border = c(0,6,1,2,3)
shape_fill = c(15,25,21,17,3)
sel_shape = c(1,3,2,4,5)
scenario_q = scenarios[[q]]; scenario_num = q

sub_df <- result_combined %>%
  filter(scenario %in% scenario_q) %>%
  mutate(
    method_g = method,
    method = paste0(method, "_", shuffle))
# Create wide DataFrames for pval, padj, and status
wide_pval <- create_wide_df(sub_df, "pval")
wide_padj <- create_wide_df(sub_df, "padj")
wide_status <- create_wide_df(sub_df, "status")
truth <- data.frame(wide_status[,1]);  colnames(truth) <- "status"; rownames(truth) <- rownames(wide_status)
# Create COBRA data and calculate performance
DF_COBRA <- COBRAData(pval = wide_pval,
                      padj = wide_padj,
                      truth = truth)
perf <- calculate_performance(DF_COBRA, binary_truth = "status")
methods_list <- c("DESpace",
                  "FindMarkers",
                  "findMarkers",
                  "pseudo.bulk_FindMarkers",
                  "spatialLIBD")
combo <- expand.grid(method = methods_list, shuffle = shuffles, stringsAsFactors = FALSE)
methods_order <- paste(combo$method, combo$shuffle, sep = "_")
methods_names <- c("DESpace2",
                   "FindMarkers (Seurat)",
                   "findMarkers (scran)",
                   "pseudo-bulk FindMarkers",
                   "spatialLIBD")

names(colors_method) <- gsub("_BayesSpace", "", names(colors_method))
names(all_colours) <- gsub("_BayesSpace", "", names(all_colours))

grad_list <- lapply(methods_list, function(m){
  base_col <- all_colours[m]
  pal <- gradient_n_pal(c(lighten(base_col, 0.6), base_col, darken(base_col, 0.6)))
  pal(seq(0, 1, length.out = length(shuffles)))
})
names(grad_list) <- methods_list

# Flatten to a named vector corresponding to method_shuffle
grad_colors <- unlist(lapply(methods_list, function(m){
  setNames(grad_list[[m]], paste(m, shuffles, sep="_"))
}))

# Prepare for plotting
cobra_plot <- prepare_data_for_plot(
  perf,
  colorscheme = c(grad_colors, "white"),
  incloverall = FALSE,
  facetted = TRUE,
  conditionalfill = FALSE
)

# Plot ROC and FDR/TPR curves
n <- length(shuffles)
(individual_v3 <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                   pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5), limits = c(0,0.5)) +
    scale_color_manual(values = grad_colors,
                       name = "",
                       breaks=methods_order,
                       labels=rep(methods_names, length(shuffles)))+
    geom_point(size = 5, aes(colour = method, shape = method
    ),
    shape = rep(c(rep(0,n), rep(c(1,6),n),
                  rep(2,n), rep(3,n)), 3), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
    scale_fill_manual(values = rep(grad_colors,3),
                      name = "",
                      breaks= methods_order,
                      labels=rep(methods_names, length(shuffles)))+
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(c(rep(15,n), rep(c(21,25),n),
                  rep(17,n), rep(3,n)), 3), 
    #rep(rep(shape_border[sel_shape],3), each = length(shuffles)), 
    stroke = 2, alpha = 0.25)+
    guides(colour = guide_legend(ncol = 1, byrow = FALSE,
                                 override.aes = list(shape = rep(shape_fill, length(shuffles)),
                                                     fill = grad_colors) ) ) +
    fig_theme + theme(legend.position = "none"
    ) )


legend1 <- readRDS(paste0(path_fig, "legend_methods.rds"))
legend2 <- readRDS(paste0(path_fig, "legend_mis_clustering_re-assigned.rds"))
max_width <- unit.pmax(legend1$widths, legend2$widths)
legend1$widths <- max_width
legend2$widths <- max_width

AA <- ggpubr::ggarrange(
  gg_fdr_ARTISTA_v3 +
    scale_y_continuous(trans = revsqrt_trans, limits = c(0.7, 1)) +
    labs(title = "Global test") +
    theme(plot.title = element_text(face = "bold")),
  
  individual_v3 +
    scale_y_continuous(trans = revsqrt_trans, limits = c(0.7, 1)) +
    labs(title = "Individual cluster") +
    theme(plot.title = element_text(face = "bold")),
  
  gridExtra::grid.arrange(legend1, legend2, ncol = 1),
  widths = c(3, 3, 1),
  ncol = 3
)

ggsave(filename = paste0(path_fig, 'Fig4.pdf'),
       plot = AA,
       device = "pdf",
       width = 16,
       height = 6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)