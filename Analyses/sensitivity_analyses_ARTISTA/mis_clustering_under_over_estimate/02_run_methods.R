rm(list = ls())
source("./sensitive_analyses/01_competitors_wrapper.R", echo=TRUE)
library(limma)
library(BiocParallel)
library(edgeR)
library(scuttle)
library(DESpace) # from 2.0
library(SingleCellExperiment)
library(edgeR)
# library(muscat)
`%notin%` <- Negate(`%in%`)

cluster_cols <- paste0("nclust_", c(3,4))
output_path <- "./Simulated_data/mis_cluster_v2/"


for (cluster_col in cluster_cols) {
  message(cluster_col)
  path <- paste0(output_path, cluster_col, "/")
  spe <- readRDS(file.path(output_path, cluster_col, "combined_simulated_Banksy_matched.rds"))
  col <- "Banksy_cluster_matched"
  ## DESpace
  res_edgeR <- dsp_test(spe,
                        cluster_col = col,
                        sample_col = "Sample",
                        condition_col = "Condition",
                        filter_gene = FALSE,
                        filter_cluster = FALSE,
                        verbose = FALSE)
  save(res_edgeR, file = file.path(output_path, cluster_col, 
                                   paste0("DESpace_global_", cluster_col, ".rda")))
  
  results_res1 <- individual_dsp(spe,
                                 cluster_col = col,
                                 sample_col = "Sample",
                                 condition_col = "Condition",
                                 filter_gene = FALSE,
                                 filter_cluster = TRUE)
  
  save(results_res1, file = file.path(output_path, cluster_col,
                                      paste0("DESpace_individual_", cluster_col, ".rda")))
  
  ## Competitors
  spe$Condition <- ifelse(spe$Condition == "Condition_1", "Condition1", "Condition2")
  run_seurat_FindMarkers(spe,
                         output_path = path,
                         cluster_method = cluster_col,
                         cluster_col = col)
  
  run_seurat_PseudoBulk_FindMarkers(spe,
                                    output_path = path,
                                    cluster_method = cluster_col,
                                    cluster_col = col)
  run_spatialLIBD_pairwise(spe,
                           output_path = path,
                           cluster_method = cluster_col,
                           cluster_col = col)
  
  run_scran_findMarkers(spe,
                        output_path = path,
                        cluster_method = cluster_col,
                        cluster_col = col)
  
}