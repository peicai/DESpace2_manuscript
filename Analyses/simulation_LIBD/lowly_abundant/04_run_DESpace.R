rm(list = ls())
load("./simulated_BayesSpace_matched.rda")
library(SingleCellExperiment)
library(dplyr)
library(edgeR)
library(muscat)
library(Seurat)
library(DESpace)
`%notin%` <- Negate(`%in%`)

sce_list_rematched <- lapply(sce_list_matched, function(x){
  sample_id <- unique(x$Sample)
  colnames(rowData(x)) <- paste0(colnames(rowData(x)), "_", sample_id)
  x
})

sce_combined = do.call(cbind, sce_list_rematched)
rm(sce_list_matched); rm(sce_list_rematched)
sce_combined_filtered <- sce_combined
level_mapping <- c("1" = "WM", "2" = "Layer6", 
                   "3" = "Layer5", "4" = "Layer4",  
                   "5" = "Layer3")

# Convert the factor levels according to the mapping
sce_combined_filtered$Cluster <-  factor(level_mapping[as.character(sce_combined_filtered$matched_cluster)])

res_edgeR <- dsp_test(spe = sce_combined_filtered,
                          cluster_col = "Cluster",
                          sample_col = "Sample",
                          condition_col = "Condition",
                          min_pct_cells = 0.5,
                          filter_gene = FALSE,
                          filter_cluster = TRUE,
                          verbose = FALSE)
save(res_edgeR, file = paste0("./DESpace_global_BayesSpace.rda"))

######################### individual test #############################
results_res1 <- individual_dsp(sce_combined_filtered,
                               cluster_col = "Cluster",
                               sample_col = "Sample",
                               condition_col = "Condition",
                               filter_gene = FALSE,
                               filter_cluster = TRUE)

save(results_res1, file = paste0(path, "DESpace_individual_BayesSpace.rda"))

#############################################################################################################
##################################################### GT #################################################
##################################################################################################################

rm(list = ls())
load("./simulated_BayesSpace_matched.rda")
library(SingleCellExperiment)
library(dplyr)
library(edgeR)
library(muscat)
library(Seurat)
library(DESpace)
`%notin%` <- Negate(`%in%`)

sce_list_rematched <- lapply(sce_list_matched, function(x){
  sample_id <- unique(x$Sample)
  colnames(rowData(x)) <- paste0(colnames(rowData(x)), "_", sample_id)
  x
})

sce_combined = do.call(cbind, sce_list_rematched)
rm(sce_list_matched); rm(sce_list_rematched)
table(sce_combined$layer_guess_reordered_droplevel)
res_edgeR <- dsp_test(spe = sce_combined,
                          cluster_col = "layer_guess_reordered_droplevel",
                          sample_col = "Sample",
                          condition_col = "Condition",
                          min_pct_cells = 0.5,
                          filter_gene = FALSE,
                          filter_cluster = TRUE,
                          verbose = FALSE)

save(res_edgeR, file = paste0("./DESpace_global_GT.rda"))

results_res1 <- individual_dsp(sce_combined,
                                cluster_col = "layer_guess_reordered_droplevel",
                                sample_col = "Sample",
                                condition_col = "Condition",
                                filter_gene = FALSE,
                                filter_cluster = TRUE)

save(results_res1, file = paste0("./DESpace_individual_GT.rda"))
