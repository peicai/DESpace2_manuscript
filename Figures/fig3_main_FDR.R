rm(list = ls())
source("./figure_theme.R", echo=TRUE)
revsqrt_trans <- scales::trans_new(
  name = "revsqrt",
  transform = function(x) -sqrt(1 - x),
  inverse = function(x) 1 - x^2
)
shape_border = c(0,6,1,2,3)
# points, fill:
shape_fill = c(15,25,21,17,3)
sel_shape = c(1,3,2,4,5)
my_theme = theme(legend.position = "none", aspect.ratio = 1,
                 plot.title = element_text(hjust = 0, face = "bold", size=13),
                 plot.subtitle = element_text(hjust = 0.5, face = "bold", size=12),
                 axis.text.x = element_text(size = 11),  
                 axis.text.y = element_text(size = 11),
                 axis.title.x = element_text(size = 12) ,
                 axis.title.y = element_text(size = 12)
)
fig_theme <- theme(strip.background = element_blank(),
                   strip.text = element_blank(),
                   axis.text.x = element_text(angle = 25, hjust=0.5, size = rel(1.5)),
                   axis.text.y=element_text(size=rel(1.5)),
                   axis.title.y = element_text(size=rel(1.5)),
                   axis.title.x = element_text(size=rel(1.5)),
                   legend.title=element_text(size=rel(1)),
                   legend.text=element_text(size=rel(1.5)),
                   legend.key.width=unit(1, "cm"),
                   legend.key = element_rect(fill = "grey90", colour = "white"),
                   legend.background = element_blank(),
                   aspect.ratio = 1, legend.position="bottom",
                   legend.box="vertical", legend.margin=margin())
########################################################################
########################## LIBD global test ############################
########################################################################
path_data <- "./data/LIBD_main/"
######################### Load data #############################
pval <- read.csv(paste0(path_data,"pval_all.csv"), row.names = 1)
padj <- read.csv(paste0(path_data,"padj_all.csv"), row.names = 1)
gene_lists <- readRDS(paste0(path_data, "geneGT.rds"))

scenarios <- list(
  list(
    gene1 = unlist(gene_lists[c(1)]),  # DSP1
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  ),
  list(
    gene1 = unlist(gene_lists[c(2)]),  # DSP2
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  ),
  list(
    gene1 = unlist(gene_lists[c(1:2)]),  # DSP1 + DSP2
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  )
)

############################## Figures ################################
gg_fdr_LIBD <- list()
q <- 3
scenario = scenarios[[q]]; scenario_num = q
gene <- c(scenario$gene1, scenario$gene2)
sub_pval <- pval[rownames(pval) %in% gene, ]
sub_padj <- padj[rownames(padj) %in% gene, ]
truth <- data.frame(status = c(rep(1, length(scenario$gene1)), rep(0, length(scenario$gene2))))
rownames(truth) <- gene
methods_order <- c("DESpace_BayesSpace",
                   "FindMarkers_BayesSpace", 
                   "findMarkers_BayesSpace", 
                   "pseudo.bulk_FindMarkers_BayesSpace",
                   "spatialLIBD_BayesSpace")
methods_names <- c("DESpace2",
                   "FindMarkers (Seurat)",
                   "findMarkers (scran)", 
                   "pseudo-bulk FindMarkers",
                   "spatialLIBD")
# Create COBRA data and calculate performance
pval <- pval[, colnames(pval) %in% methods_order]
padj <- padj[, colnames(padj) %in% methods_order]

DF_COBRA <- COBRAData(pval = data.frame(pval), 
                      padj = data.frame(padj), 
                      truth = truth)
perf <- calculate_performance(DF_COBRA, binary_truth = "status")

# Prepare for plotting
cobra_plot <- prepare_data_for_plot(
  perf,
  colorscheme = all_colours,
  incloverall = FALSE,
  facetted = TRUE,
  conditionalfill = FALSE
)

# Plot ROC and FDR/TPR curves
all_colours <- all_colours[methods_order]
(gg_fdr_LIBD_1 <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                   pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.4), limits = c(0,0.4)) +
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_order,
                       labels=methods_names)+
    guides(colour = guide_legend(ncol = 3, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = all_colours) ) ) +
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],3), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
    scale_fill_manual(values = rep(all_colours,3),
                      name = "",
                      breaks= methods_order,
                      labels= methods_names)+
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_fill[sel_shape],3), 
    stroke = 2, alpha = 0.25)+
    fig_theme )


########################################################################
########################## LIBD individual test ############################
########################################################################
path_data <- "./data/LIBD_individual/"
######################### Load data #############################
result_combined <- read.csv(paste0(path_data,"results.csv"), row.names = 1)
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
############################## Figures ################################
scenario_q <- c("DSP1", "DSP2", "NULL1", "NULL2")
sub_df <- result_combined %>% filter(scenario %in% scenario_q) 
sub_df_list <- lapply(c("pval", "padj", "status"), function(col) {
  df <- data.frame(sub_df[[col]], row.names = rownames(sub_df))
  colnames(df) <- col
  return(df)
})

# Create wide DataFrames for pval, padj, and status
wide_pval <- create_wide_df(sub_df, "pval")
wide_padj <- create_wide_df(sub_df, "padj")
wide_status <- create_wide_df(sub_df, "status") 
truth <- data.frame(wide_status[,1]);  colnames(truth) <- "status"; rownames(truth) <- rownames(wide_status)
# Create COBRA data and calculate performance
wide_pval <- wide_pval[, colnames(wide_pval) %in% methods_order]
wide_padj <- wide_padj[, colnames(wide_padj) %in% methods_order]

DF_COBRA <- COBRAData(pval = wide_pval, 
                      padj = wide_padj, 
                      truth = truth)
perf <- calculate_performance(DF_COBRA, binary_truth = "status")

# Prepare for plotting
cobra_plot <- prepare_data_for_plot(
  perf,
  colorscheme = all_colours,
  incloverall = FALSE,
  facetted = TRUE,
  conditionalfill = FALSE
)


(gg_fdr_LIBD_2 <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                  pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.4), limits = c(0,0.4)) +
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_order,
                       labels=methods_names
    )+
    guides(colour = guide_legend(ncol = 3, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = all_colours) ) ) +
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],3), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
    scale_fill_manual(values = rep(all_colours,3),
                      name = "",
                      breaks= methods_order,
                      labels= methods_names
    )+
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_fill[sel_shape],3), 
    stroke = 2, alpha = 0.25)+
    fig_theme )


##################### combine TPR vs. FDR curve #################
legend <- ggpubr::get_legend(gg_fdr_LIBD_1 +
                               guides(colour = guide_legend(nrow = 2, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = all_colours) ) ) +
                               theme(legend.key.height=unit(1.5,"line"),
                                     legend.text=element_text(size=13),
                                     legend.key.width=unit(2,"line")))

subtitles <- c("global test", "individual cluster")
# Apply transformation and formatting in a loop  
(FDR_LIBD = ggpubr::ggarrange( gg_fdr_LIBD_1 + scale_y_continuous(trans = revsqrt_trans, limits = c(0, 1)) +  
                                 labs(title = "LIBD", subtitle = subtitles[1]) + my_theme,
                               gg_fdr_LIBD_2 + scale_y_continuous(trans = revsqrt_trans, limits = c(0, 1)) +  
                                 labs(title = "", subtitle = subtitles[2]) + my_theme
                               
))


########################################################################
########################## ARTISTA global test ############################
########################################################################
path_data <- "./data/ARTISTA_main/"
######################### Load data #############################
pval <- read.csv(paste0(path_data,"pval_all.csv"), row.names = 1)
padj <- read.csv(paste0(path_data,"padj_all.csv"), row.names = 1)
gene_lists <- readRDS(paste0(path_data, "geneGT.rds"))

scenarios <- list(
  list(
    gene1 = unlist(gene_lists[c(1)]),  # DSP1
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  ),
  list(
    gene1 = unlist(gene_lists[c(2)]),  # DSP2
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  ),
  list(
    gene1 = unlist(gene_lists[c(1:2)]),  # DSP1 + DSP2
    gene2 = unlist(gene_lists[c(3:4)])  # NULL
  )
)
names(colors_method) <- gsub("_BayesSpace", "_Banksy", names(colors_method))
names(all_colours) <- gsub("_BayesSpace", "_Banksy", names(all_colours))
############################## Figures ################################
gg_fdr_ARTISTA <- list()
q <- 3
scenario = scenarios[[q]]; scenario_num = q
gene <- c(scenario$gene1, scenario$gene2)
sub_pval <- pval[rownames(pval) %in% gene, ]
sub_padj <- padj[rownames(padj) %in% gene, ]
truth <- data.frame(status = c(rep(1, length(scenario$gene1)), rep(0, length(scenario$gene2))))
rownames(truth) <- gene
methods_order <- c("DESpace_Banksy",
                   "FindMarkers_Banksy",
                   "findMarkers_Banksy", 
                   "pseudo.bulk_FindMarkers_Banksy",
                   "spatialLIBD_Banksy")
methods_names <- c("DESpace2",
                   "FindMarkers (Seurat)",
                   "findMarkers (scran)", 
                   "pseudo-bulk FindMarkers",
                   "spatialLIBD")
# Create COBRA data and calculate performance
pval <- pval[, colnames(pval) %in% methods_order]
padj <- padj[, colnames(padj) %in% methods_order]

DF_COBRA <- COBRAData(pval = data.frame(pval), 
                      padj = data.frame(padj), 
                      truth = truth)
perf <- calculate_performance(DF_COBRA, binary_truth = "status")

# Prepare for plotting
cobra_plot <- prepare_data_for_plot(
  perf,
  colorscheme = all_colours,
  incloverall = FALSE,
  facetted = TRUE,
  conditionalfill = FALSE
)

# Plot ROC and FDR/TPR curves
(gg_fdr_ARTISTA_1 <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                      pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.4), limits = c(0,0.4)) +
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_order,
                       labels=methods_names)+
    guides(colour = guide_legend(ncol = 3, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = all_colours) ) ) +
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],3), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
    scale_fill_manual(values = rep(all_colours,3),
                      name = "",
                      breaks= methods_order,
                      labels= methods_names)+
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_fill[sel_shape],3), 
    stroke = 2, alpha = 0.25)+
    fig_theme + theme(legend.position = "right", legend.direction = "horizontal") )

########################## ARTISTA individual test ############################
########################################################################
path_data <- "./data/ARTISTA_individual/"
######################### Load data #############################
result_combined <- read.csv(paste0(path_data,"results.csv"), row.names = 1)
result_combined <- result_combined %>% filter(method != "DESpace_GT")
result_combined$method <- gsub("_Banksy", "", result_combined$method)
scenarios <- list(
  c("DSP1", "NULL1", "NULL2"),
  c("DSP2", "NULL1", "NULL2"),
  c("DSP1", "DSP2", "NULL1", "NULL2")
)
names(colors_method) <- gsub("_BayesSpace", "", names(colors_method))
names(all_colours) <- gsub("_BayesSpace", "", names(all_colours))
methods_order <- c("DESpace",
                   "FindMarkers", 
                   "findMarkers", 
                   "pseudo.bulk_FindMarkers",
                   "spatialLIBD")

############################## Figures ################################
q <- 3
scenario_q = scenarios[[q]]; scenario_num = q
sub_df <- result_combined %>% filter(scenario %in% scenario_q) 
# Create wide DataFrames for pval, padj, and status
wide_pval <- create_wide_df(sub_df, "pval")
wide_padj <- create_wide_df(sub_df, "padj")
wide_status <- create_wide_df(sub_df, "status") 
truth <- data.frame(wide_status[,1]);  colnames(truth) <- "status"; rownames(truth) <- rownames(wide_status)
# Create COBRA data and calculate performance
wide_pval <- wide_pval[,colnames(wide_pval) %in% methods_order]
wide_padj <- wide_padj[,colnames(wide_padj) %in% methods_order]

DF_COBRA <- COBRAData(pval = wide_pval, 
                      padj = wide_padj, 
                      truth = truth)
perf <- calculate_performance(DF_COBRA, binary_truth = "status")

# Prepare for plotting
cobra_plot <- prepare_data_for_plot(
  perf,
  colorscheme = all_colours,
  incloverall = FALSE,
  facetted = TRUE,
  conditionalfill = FALSE
)

# Plot ROC and FDR/TPR curves
(gg_fdr_ARTISTA_2 <- plot_fdrtprcurve(cobra_plot,plottype = c("points"),
                                      pointsize = 0, linewidth = 2)+
    scale_x_sqrt( breaks = c(0.00, 0.01, 0.05, 0.1, 0.2, 0.4), limits = c(0,0.4)) +
    scale_color_manual(values = all_colours,
                       name = "",
                       breaks=methods_order,
                       labels=methods_names
    )+
    guides(colour = guide_legend(ncol = 3, byrow = FALSE,
                                 override.aes = list(shape = shape_fill,
                                                     fill = all_colours) ) ) +
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_border[sel_shape],3), 
    stroke = 2, 
    alpha = 1) + # stroke = line width
    scale_fill_manual(values = rep(all_colours,3),
                      name = "",
                      breaks= methods_order,
                      labels= methods_names
    )+
    geom_point(size = 5, aes(fill = method, colour = method, shape = method
    ),
    shape = rep(shape_fill[sel_shape],3), 
    stroke = 2, alpha = 0.25)+
    fig_theme )

##################### combine TPR vs. FDR curve #################
legend <- ggpubr::get_legend(gg_fdr_ARTISTA_1 +
                               guides(colour = guide_legend(nrow = 10, byrow = FALSE,
                                                            override.aes = list(shape = shape_fill,
                                                                                fill = all_colours) ) ) +
                               theme(legend.key.height=unit(1.5,"line"),
                                     legend.text=element_text(size=11),
                                     legend.key.width=unit(2,"line")))



(FDR_ARTISTA = ggpubr::ggarrange( gg_fdr_ARTISTA_1 + scale_y_continuous(trans = revsqrt_trans,
                                                                        limits = c(0, 1)) +  
                                    labs(title = "", subtitle = "ARTISTA") + my_theme + 
                                    theme(
                                      plot.subtitle = element_text(hjust = 0, face = "bold", size=13),
                                      plot.margin = margin(5, 5, 5, 5)
                                    ),
                                  gg_fdr_ARTISTA_2 + scale_y_continuous(trans = revsqrt_trans,
                                                                        limits = c(0, 1)) +  
                                    labs(title = "", subtitle = "") + my_theme 
))

(AA = ggpubr::ggarrange(FDR_LIBD, FDR_ARTISTA,
                        heights=c(1,1),
                        ncol = 1 ))
(AAA = ggpubr::ggarrange(AA, legend,
                         widths =c(3,1),
                         ncol = 2 ))
path_fig <- "./figure/"

ggsave(filename = paste0(path_fig, 'main_FDR.pdf'),
       plot = AAA,
       device = "pdf",
       width = 12,
       height = 8,
       units = "in",
       dpi = 300,
       limitsize = TRUE)

