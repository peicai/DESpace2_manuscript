rm(list = ls())
load("./simulated_BayesSpace_matched.rda")
library(SingleCellExperiment)
library(dplyr)
library(edgeR)
library(muscat)
`%notin%` <- Negate(`%in%`)
sce_list_rematched <- lapply(sce_list_matched, function(x){
  sample_id <- unique(x$Sample)
  colnames(rowData(x)) <- paste0(colnames(rowData(x)), "_", sample_id)
  x
})

sce_combined = do.call(cbind, sce_list_rematched)
rm(sce_list_matched); rm(sce_list_rematched)
source("./01_competitors_wrapper.R", echo=TRUE)
sce_combined_filtered <- sce_combined
level_mapping <- c("1" = "WM", "2" = "Layer6", 
                   "3" = "Layer5", "4" = "Layer4",  
                   "5" = "Layer3")

# Convert the factor levels according to the mapping
sce_combined_filtered$Cluster <-  factor(level_mapping[as.character(sce_combined_filtered$matched_cluster)])


run_seurat_FindMarkers(sce_combined_filtered,
                          output_path = "./",
                          cluster_method = "BayesSpace",
                          cluster_col = "Cluster")
  
run_seurat_PseudoBulk_FindMarkers(sce_combined_filtered,
                          output_path = "./",
                          cluster_method = "BayesSpace",
                          cluster_col = "Cluster")

run_LIBD_registration1(sce_combined_filtered,
                      output_path = "./",
                      cluster_method = "BayesSpace",
                      cluster_col = "Cluster")

run_scran_findMarkers(sce_combined_filtered,
                      output_path = "./",
                      cluster_method = "BayesSpace",
                      cluster_col = "Cluster")

########################### GT #############################
sce_combined_filtered <- sce_combined
table(sce_combined_filtered$layer_guess_reordered_droplevel)
run_seurat_FindMarkers(sce_combined_filtered,
                          output_path = "./",
                          cluster_method = "GT",
                          cluster_col = "layer_guess_reordered_droplevel")

run_seurat_PseudoBulk_FindMarkers(sce_combined_filtered,
                                  output_path = "./",
                                  cluster_method = "GT",
                                  cluster_col = "layer_guess_reordered_droplevel")

run_LIBD_registration1(sce_combined_filtered,
                       output_path = "./",
                       cluster_method = "GT",
                       cluster_col = "layer_guess_reordered_droplevel")

run_scran_findMarkers(sce_combined_filtered,
                      output_path = "./",
                      cluster_method = "GT",
                      cluster_col = "layer_guess_reordered_droplevel")

load("./seurat_PseudoBulk_FindMarkers_GT.rda")
all_genes <- rownames(sce_combined)
miss_genes <- all_genes[all_genes %notin% rownames(results_res$Layer4)]
dummy_rows <- data.frame(
  do.call(rbind, lapply(miss_genes, function(gene) c(p_val = 1, 
                                                     avg_log2FC = 1, 
                                                     pct.1 = 1, 
                                                     pct.2 = 1, 
                                                     p_val_adj = 1))),
  stringsAsFactors = FALSE
)

rownames(dummy_rows) <- miss_genes
results_res$Layer4 <- rbind(results_res$Layer4, dummy_rows)

miss_genes <- all_genes[all_genes %notin% rownames(results_res$WM)]
dummy_rows <- data.frame(
  do.call(rbind, lapply(miss_genes, function(gene) c(p_val = 1, 
                                                     avg_log2FC = 1, 
                                                     pct.1 = 1, 
                                                     pct.2 = 1, 
                                                     p_val_adj = 1))),
  stringsAsFactors = FALSE
)

rownames(dummy_rows) <- miss_genes
results_res$WM <- rbind(results_res$WM, dummy_rows)
save(results_res, file ="./seurat_PseudoBulk_FindMarkers_GT.rda")