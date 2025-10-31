rm(list = ls())
source("./figure_theme.R", echo=TRUE)
path_fig = "./figure/"
path_data <- "./data/"
shape_border = c(0,6,1,2,3)
shape_fill = c(15,25,21,17,3)
sel_shape = c(1,3,2,4,5)
fig_theme <- theme(strip.background = element_blank(),
                   strip.text = element_blank(),
                   axis.text.x = element_text(angle = 25, hjust=0.5, size = rel(1)),
                   axis.text.y=element_text(size=rel(1)),
                   axis.title.y = element_text(size=rel(1)),
                   axis.title.x = element_text(size=rel(1)),
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
    gene1 = unlist(gene_lists[c(1:2)]),  # DSP1 + DSP2
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  )
names(colors_method) <- gsub("_BayesSpace", "", names(colors_method))
names(all_colours) <- gsub("_BayesSpace", "", names(all_colours))
############################## Figures ################################
library(dplyr)
shuffles <- c("0", "001", "005", "01", "02")
padj_list <- list()
pval_list <- list()
for (s in shuffles) {
  cols_s <- grep(paste0("_", s, "$"), colnames(padj), value = TRUE)
  df <- padj[, cols_s, drop = FALSE]
  colnames(df) <- gsub(paste0("_", s, "$"), "", colnames(df))
  padj_list[[s]] <- df
  
  cols_s <- grep(paste0("_", s, "$"), colnames(pval), value = TRUE)
  df <- pval[, cols_s, drop = FALSE]
  colnames(df) <- gsub(paste0("_", s, "$"), "", colnames(df))
  pval_list[[s]] <- df
}

gg_fdr_ARTISTA <- list()
gene <- c(scenario$gene1, scenario$gene2)
all_colours[["spatialLIBD_pairwise"]] <- all_colours[["spatialLIBD"]]
title_shuffle <- c("0%", "1%", "5%", "10%", "20%")
for(i in seq(shuffles)){
  s <- shuffles[i]
  sub_pval <- pval_list[[s]][rownames(pval_list[[s]]) %in% gene, ]
  sub_padj <- padj_list[[s]][rownames(padj_list[[s]]) %in% gene, ]
  truth <- data.frame(status = c(rep(1, length(scenario$gene1)), rep(0, length(scenario$gene2))))
  rownames(truth) <- gene
  methods_order <- c("DESpace",
                     "FindMarkers", 
                     "findMarkers",
                     "pseudo.bulk_FindMarkers",
                     "spatialLIBD_pairwise")
  methods_names <- c("DESpace2",
                     "FindMarkers (Seurat)",
                     "findMarkers (scran)",
                     "pseudo-bulk FindMarkers",
                     "spatialLIBD")
  all_colours_update <- all_colours[methods_order]
  DF_COBRA <- COBRAData(pval = data.frame(sub_pval), 
                        padj = data.frame(sub_padj), 
                        truth = truth)
  perf <- calculate_performance(DF_COBRA, binary_truth = "status")
  cobra_plot <- prepare_data_for_plot(
    perf,
    colorscheme = c(all_colours_update, "white"),
    incloverall = FALSE,
    facetted = TRUE,
    conditionalfill = FALSE
  )
  
  # Plot ROC and FDR/TPR curves
  (gg_fdr_ARTISTA[[s]] <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                           pointsize = 0, linewidth = 2)+
      scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.4), limits = c(0,0.4)) +
      scale_color_manual(values = all_colours,
                         name = "",
                         breaks=methods_order,
                         labels=methods_names)+
      guides(colour = guide_legend(ncol = 5, byrow = FALSE,
                                   override.aes = list(shape = shape_fill,
                                                       fill = all_colours_update) ) ) +
      geom_point(size = 5, aes(fill = method, colour = method, shape = method
      ),
      shape = rep(shape_border[sel_shape],3), 
      stroke = 2, 
      alpha = 1) + 
      scale_fill_manual(values = rep(all_colours_update,3),
                        name = "",
                        breaks= methods_order,
                        labels= methods_names)+
      geom_point(size = 5, aes(fill = method, colour = method, shape = method
      ),
      shape = rep(shape_fill[sel_shape],3), 
      stroke = 2, alpha = 0.25)+
      fig_theme + theme(legend.position = "none", legend.direction = "vertical"
      ) +
      scale_y_continuous(trans = revsqrt_trans,
                         limits = c(0, 1)) +
      labs(title = title_shuffle[i])
  ) 
}

FDR_ARTISTA_global <- do.call(ggarrange, c(gg_fdr_ARTISTA, ncol = 5, nrow = 1))

########################################################################
########################## ARTISTA individual test ############################
########################################################################
######################### Load data #############################
library(readr)
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
library(dplyr)
shuffles <- c("0", "001", "005", "01", "02")
q <- 3
source("~/Desktop/PhD_Project/DESpace_extension/manuscript_code/scripts/figure_theme.R", echo=TRUE)
shape_border = c(0,6,1,2,3)
shape_fill = c(15,25,21,17,3)
sel_shape = c(1,3,2,4,5)
scenario_q = scenarios[[q]]; scenario_num = q
############################## Figures ################################
names(colors_method) <- gsub("_BayesSpace", "", names(colors_method))
names(all_colours) <- gsub("_BayesSpace", "", names(all_colours))

gg_fdr_ARTISTA2 <- list()
for(i in seq(shuffles)){
  s <- shuffles[i]
  sub_df <- result_combined %>% filter(scenario %in% scenario_q,
                                       shuffle == s) 
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
  methods_order <- c("DESpace",
                     "FindMarkers", 
                     "findMarkers",
                     "pseudo.bulk_FindMarkers",
                     "spatialLIBD")
  methods_names <- c("DESpace2",
                     "FindMarkers (Seurat)",
                     "findMarkers (scran)",
                     "pseudo-bulk FindMarkers",
                     "spatialLIBD")
  all_colours_update <- all_colours[methods_order]
  # Prepare for plotting
  cobra_plot <- prepare_data_for_plot(
    perf,
    colorscheme = c(all_colours_update, "white"),
    incloverall = FALSE,
    facetted = TRUE,
    conditionalfill = FALSE
  )
  (gg_fdr_ARTISTA2[[s]] <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                            pointsize = 0, linewidth = 2)+
      scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.5), limits = c(0,0.5)) +
      scale_color_manual(values = all_colours,
                         name = "",
                         breaks=methods_order,
                         labels=methods_names)+
      guides(colour = guide_legend(ncol = 1, byrow = FALSE,
                                   override.aes = list(shape = shape_fill,
                                                       fill = all_colours_update) ) ) +
      geom_point(size = 5, aes(fill = method, colour = method, shape = method
      ),
      shape = rep(shape_border[sel_shape],3), 
      stroke = 2, 
      alpha = 1) + # stroke = line width
      scale_fill_manual(values = rep(all_colours_update,3),
                        name = "",
                        breaks= methods_order,
                        labels= methods_names)+
      geom_point(size = 5, aes(fill = method, colour = method, shape = method
      ),
      shape = rep(shape_fill[sel_shape],3), 
      stroke = 2, alpha = 0.25)+
      fig_theme + theme(legend.position = "none"
      )+
      scale_y_continuous(trans = revsqrt_trans,
                         limits = c(0, 1)) +
      labs(title = title_shuffle[i])
  )
  
}

FDR_ARTISTA_individual <- do.call(ggarrange, c(gg_fdr_ARTISTA2, ncol = 5, nrow = 1))

library(ggpubr)

FDR_ARTISTA_individual2 <- annotate_figure(
  FDR_ARTISTA_individual,
  top =  text_grob(
    "Individual cluster",
    size = 14,
    face = "bold",  
    hjust = 0,       
    x = 0.04,         
    y = 0.4     
  )
)

FDR_ARTISTA_global2 <- annotate_figure(
  FDR_ARTISTA_global,
  top =  text_grob(
    "Global test",
    size = 14,
    face = "bold",  
    hjust = 0,       
    x = 0.038,        
    y = 0.4       
  )
)
source("./figure_theme.R", echo=TRUE)
legend <- readRDS(paste0(path_fig, "legend_methods_bottom.rds"))
(AA = ggpubr::ggarrange(FDR_ARTISTA_global2, 
                        FDR_ARTISTA_individual2,
                        legend,
                        heights = c(5,5, 1),
                        ncol = 1 ))

ggsave(filename = paste0(path_fig, 'supp10.pdf'),
       plot = AA,
       device = "pdf",
       width = 15,
       height = 7.6,
       units = "in",
       dpi = 300,
       limitsize = TRUE)